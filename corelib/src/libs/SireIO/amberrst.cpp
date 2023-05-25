/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2023  Christopher Woods
  *
  *  This program is free software; you can redistribute it and/or modify
  *  it under the terms of the GNU General Public License as published by
  *  the Free Software Foundation; either version 3 of the License, or
  *  (at your option) any later version.
  *
  *  This program is distributed in the hope that it will be useful,
  *  but WITHOUT ANY WARRANTY; without even the implied warranty of
  *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *  GNU General Public License for more details.
  *
  *  You should have received a copy of the GNU General Public License
  *  along with this program; if not, write to the Free Software
  *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
  *
  *  For full details of the license please see the COPYING file
  *  that should have come with this distribution.
  *
  *  You can contact the authors via the website
  *  at http://sire.openbiosim.org
  *
\*********************************************/

#include "SireIO/amberrst.h"
#include "SireIO/netcdffile.h"
#include "SireIO/xdrfile.h"

#include "SireSystem/system.h"

#include "SireMol/atomcoords.h"
#include "SireMol/atomvelocities.h"
#include "SireMol/core.h"
#include "SireMol/mgname.h"
#include "SireMol/molecule.h"
#include "SireMol/moleditor.h"
#include "SireMol/molidx.h"
#include "SireMol/trajectory.h"

#include "SireVol/periodicbox.h"
#include "SireVol/triclinicbox.h"

#include "SireUnits/dimensions.h"
#include "SireUnits/units.h"

#include "SireBase/generalunitproperty.h"
#include "SireBase/parallel.h"
#include "SireBase/stringproperty.h"
#include "SireBase/timeproperty.h"
#include "SireBase/progressbar.h"
#include "SireBase/releasegil.h"

#include "SireIO/errors.h"
#include "SireError/errors.h"

#include "SireStream/shareddatastream.h"

#include "tostring.h"

#include <QDebug>

using namespace SireIO;
using namespace SireIO::detail;
using namespace SireMaths;
using namespace SireMol;
using namespace SireBase;
using namespace SireSystem;
using namespace SireVol;
using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireStream;

//////////
////////// Implementation of AmberRstFrameBuffer
//////////

class AmberRstFrameBuffer
{
public:
    AmberRstFrameBuffer(const Frame &frame)
        : current_frame(0), natoms(0)
    {
        if (frame.nAtoms() == 0)
            return;

        natoms = frame.nAtoms();

        if (frame.hasCoordinates())
            coords = QVector<double>(3 * natoms);

        if (frame.hasVelocities())
            vels = QVector<double>(3 * natoms);

        if (frame.hasForces())
            frcs = QVector<double>(3 * natoms);

        if (frame.space().isPeriodic())
            spce = QVector<double>(6);
    }

    ~AmberRstFrameBuffer()
    {
    }

    int current_frame;
    int natoms;

    QVector<double> coords;
    QVector<double> vels;
    QVector<double> frcs;
    QVector<double> box;
    QVector<double> spce;

    double time;
};

//////////
////////// Implementation of AmberRstFile
//////////

class AmberRstFile : public NetCDFFile
{
public:
    AmberRstFile();
    AmberRstFile(const QString &filename);
    ~AmberRstFile();

    bool open(QIODevice::OpenMode mode = QIODevice::ReadOnly,
              bool use_64bit_offset = true, bool use_netcdf4 = true);

    SireMol::Frame readFrame(int i, bool use_parallel = true) const;
    void writeFrame(const SireMol::Frame &frame,
                    bool use_parallel = true);

    int nAtoms() const;
    int nFrames() const;

private:
    void _lkr_reset();
    void _lkr_reindexFrames();
    void _lkr_readFrameIntoBuffer(int i);
    void _lkr_writeBufferToFile();

    /** The current frame that has been read into the buffer */
    std::shared_ptr<AmberRstFrameBuffer> frame_buffer;

    /** The application that created this file */
    QString application;

    /** The specific program that created this file */
    QString program;

    /** The version of the program */
    QString program_version;

    /** Title for the file */
    QString title;

    /** The number of atoms in the frame - we assume all
     *  frames have the same number of atoms
     */
    qint64 natoms;

    /** The number of frames in the file */
    qint64 nframes;

    /** Whether or not this is a restart file */
    bool is_restart;
};

AmberRstFile::AmberRstFile()
    : NetCDFFile(),
      natoms(0), nframes(0)
{
}

AmberRstFile::AmberRstFile(const QString &filename)
    : NetCDFFile(filename),
      natoms(0), nframes(0)
{
}

AmberRstFile::~AmberRstFile()
{
}

int AmberRstFile::nAtoms() const
{
    return natoms;
}

int AmberRstFile::nFrames() const
{
    return nframes;
}

void AmberRstFile::_lkr_reset()
{
    frame_buffer.reset();

    natoms = 0;
    nframes = 0;
}

void AmberRstFile::_lkr_writeBufferToFile()
{
}

void AmberRstFile::writeFrame(const Frame &frame, bool use_parallel)
{
    // create a frame buffer for this frame
    std::shared_ptr<AmberRstFrameBuffer> buffer(new AmberRstFrameBuffer(frame));

    // copy the space info into the buffer

    buffer->time = frame.time().to(picosecond);

    const int natoms = frame.nAtoms();

    if (natoms != buffer->natoms)
        throw SireError::program_bug(QObject::tr(
                                         "Memory not allocated!"),
                                     CODELOC);

    auto copy_coordinates = [&]()
    {
        if (not frame.hasCoordinates())
            return;

        const auto coords_data = frame.coordinates().constData();
        auto c_data = buffer->coords.data();

        if (c_data == 0)
            return;

        if (use_parallel)
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, natoms), [&](tbb::blocked_range<int> r)
                              {
                for (int i=r.begin(); i<r.end(); ++i)
                {
                    const auto &value = coords_data[i];

                    c_data[(3*i) + 0] = value.x();
                    c_data[(3*i) + 1] = value.y();
                    c_data[(3*i) + 2] = value.z();
                } });
        }
        else
        {
            for (int i = 0; i < natoms; ++i)
            {
                const auto &value = coords_data[i];

                c_data[(3 * i) + 0] = value.x();
                c_data[(3 * i) + 1] = value.y();
                c_data[(3 * i) + 2] = value.z();
            }
        }
    };

    auto copy_velocities = [&]()
    {
        if (not frame.hasVelocities())
            return;

        const auto vels_data = frame.velocities().constData();
        auto v_data = buffer->vels.data();

        if (v_data == 0)
            return;

        if (use_parallel)
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, natoms), [&](tbb::blocked_range<int> r)
                              {
                for (int i=r.begin(); i<r.end(); ++i)
                {
                    const auto &value = vels_data[i];

                    v_data[(3*i) + 0] = value.x().to(nanometer / picosecond);
                    v_data[(3*i) + 1] = value.y().to(nanometer / picosecond);
                    v_data[(3*i) + 2] = value.z().to(nanometer / picosecond);
                } });
        }
        else
        {
            for (int i = 0; i < natoms; ++i)
            {
                const auto &value = vels_data[i];

                v_data[(3 * i) + 0] = value.x().to(nanometer / picosecond);
                v_data[(3 * i) + 1] = value.y().to(nanometer / picosecond);
                v_data[(3 * i) + 2] = value.z().to(nanometer / picosecond);
            }
        }
    };

    auto copy_forces = [&]()
    {
        if (not frame.hasForces())
            return;

        const auto frcs_data = frame.forces().constData();
        auto f_data = buffer->frcs.data();

        if (f_data == 0)
            return;

        if (use_parallel)
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, natoms), [&](tbb::blocked_range<int> r)
                              {
                for (int i=r.begin(); i<r.end(); ++i)
                {
                    const auto &value = frcs_data[i];

                    f_data[(3*i) + 0] = value.x().to(kilojoule / nanometer);
                    f_data[(3*i) + 1] = value.y().to(kilojoule / nanometer);
                    f_data[(3*i) + 2] = value.z().to(kilojoule / nanometer);
                } });
        }
        else
        {
            for (int i = 0; i < natoms; ++i)
            {
                const auto &value = frcs_data[i];

                f_data[(3 * i) + 0] = value.x().to(kilojoule / nanometer);
                f_data[(3 * i) + 1] = value.y().to(kilojoule / nanometer);
                f_data[(3 * i) + 2] = value.z().to(kilojoule / nanometer);
            }
        }
    };

    if (use_parallel)
    {
        tbb::parallel_invoke([&]()
                             { copy_coordinates(); },
                             [&]()
                             { copy_velocities(); },
                             [&]()
                             { copy_forces(); });
    }
    else
    {
        copy_coordinates();
        copy_velocities();
        copy_forces();
    }

    QMutexLocker lkr(NetCDFFile::globalMutex());

    // replace the current buffer with this new buffer
    frame_buffer.reset();
    frame_buffer = buffer;

    // now write this to the file
    this->_lkr_writeBufferToFile();
}

void AmberRstFile::_lkr_readFrameIntoBuffer(int frame)
{
    if (frame < 0 or frame >= nframes)
    {
        throw SireError::invalid_index(QObject::tr(
                                           "Cannot read frame %1 as the number of frames is %2.")
                                           .arg(frame)
                                           .arg(nframes),
                                       CODELOC);
    }

    qDebug() << "READ FRAME" << frame;

    const auto dimensions = this->_lkr_getVariablesInfo();

    const auto coordinates = dimensions.value("coordinates");
    const auto slab = coordinates.hyperslab();

    qDebug() << coordinates.toString();
    qDebug() << slab.toString();

    auto coordinates_data = this->_lkr_read(coordinates, slab);

    qDebug() << coordinates_data.toString();
    qDebug() << Sire::toString(coordinates_data.toFloatArray());
}

Frame AmberRstFile::readFrame(int i, bool use_parallel) const
{
    AmberRstFile *nonconst_this = const_cast<AmberRstFile *>(this);

    QMutexLocker lkr(NetCDFFile::globalMutex());

    // load the frame into the buffer
    nonconst_this->_lkr_readFrameIntoBuffer(i);

    if (frame_buffer.get() == 0)
        // no frame has been loaded?
        return Frame();

    // copy it out from the buffer to local storage
    const auto natoms = frame_buffer->natoms;

    SpacePtr space;
    auto time = double(frame_buffer->time) * picosecond;

    // and also the coords / vels / forces data
    QVector<double> coords = frame_buffer->coords;
    QVector<double> vels = frame_buffer->vels;
    QVector<double> frcs = frame_buffer->frcs;

    // we've finished with the buffer - can release the mutex
    lkr.unlock();

    // now need to convert the coords, vels and frcs into the right units
    // Can do this in parallel if allowed :-)
    QVector<Vector> c;
    QVector<Velocity3D> v;
    QVector<Force3D> f;

    auto copy_coordinates = [&]()
    {
        if (coords.isEmpty())
            return;

        c = QVector<Vector>(natoms);
        auto c_data = c.data();
        auto coords_data = coords.constData();

        if (use_parallel)
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, natoms), [&](const tbb::blocked_range<int> &r)
                              {
                for (int i=r.begin(); i<r.end(); ++i)
                {
                    c_data[i] = Vector(coords_data[(3 * i) + 0],
                                       coords_data[(3 * i) + 1],
                                       coords_data[(3 * i) + 2]);
                } });
        }
        else
        {
            for (int i = 0; i < natoms; ++i)
            {
                c_data[i] = Vector(coords_data[(3 * i) + 0],
                                   coords_data[(3 * i) + 1],
                                   coords_data[(3 * i) + 2]);
            }
        }
    };

    auto copy_velocities = [&]()
    {
        if (vels.isEmpty())
            return;

        v = QVector<Velocity3D>(natoms);
        auto v_data = v.data();
        auto vels_data = vels.constData();

        const auto units_to_internal = (1 * nanometer / picosecond);

        if (use_parallel)
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, natoms), [&](const tbb::blocked_range<int> &r)
                              {
                for (int i=r.begin(); i<r.end(); ++i)
                {
                    v_data[i] = Velocity3D(vels_data[(3 * i) + 0] * units_to_internal,
                                           vels_data[(3 * i) + 1] * units_to_internal,
                                           vels_data[(3 * i) + 2] * units_to_internal);
                } });
        }
        else
        {
            for (int i = 0; i < natoms; ++i)
            {
                v_data[i] = Velocity3D(vels_data[(3 * i) + 0] * units_to_internal,
                                       vels_data[(3 * i) + 1] * units_to_internal,
                                       vels_data[(3 * i) + 2] * units_to_internal);
            }
        }
    };

    auto copy_forces = [&]()
    {
        if (frcs.isEmpty())
            return;

        f = QVector<Force3D>(natoms);
        auto f_data = f.data();
        auto frcs_data = frcs.constData();

        const auto units_to_internal = (kilojoule / nanometer);

        if (use_parallel)
        {
            tbb::parallel_for(tbb::blocked_range<int>(0, natoms), [&](const tbb::blocked_range<int> &r)
                              {
                for (int i=r.begin(); i<r.end(); ++i)
                {
                    f_data[i] = Force3D(frcs_data[(3 * i) + 0] * units_to_internal,
                                        frcs_data[(3 * i) + 1] * units_to_internal,
                                        frcs_data[(3 * i) + 2] * units_to_internal);
                } });
        }
        else
        {
            for (int i = 0; i < natoms; ++i)
            {
                f_data[i] = Force3D(frcs_data[(3 * i) + 0] * units_to_internal,
                                    frcs_data[(3 * i) + 1] * units_to_internal,
                                    frcs_data[(3 * i) + 2] * units_to_internal);
            }
        }
    };

    if (use_parallel)
    {
        tbb::parallel_invoke([&]()
                             { copy_coordinates(); },
                             [&]()
                             { copy_velocities(); },
                             [&]()
                             { copy_forces(); });
    }
    else
    {
        copy_coordinates();
        copy_velocities();
        copy_forces();
    }

    // NEED TO CONSTRUCT THE SPACE

    return Frame(c, v, f, Cartesian(), time);
}

bool AmberRstFile::open(QIODevice::OpenMode mode, bool use_64bit_offset, bool use_netcdf4)
{
    QMutexLocker lkr(NetCDFFile::globalMutex());

    if (not this->_lkr_open(mode, use_64bit_offset, use_netcdf4))
        return false;

    if (mode != QIODevice::ReadOnly)
    {
        // we only want to write, and so we don't need to check the file
        return true;
    }

    // now get the attributes from the file
    QString conventions;

    try
    {
        conventions = this->_lkr_getStringAttribute("Conventions");
    }
    catch (...)
    {
        this->_lkr_close();
        throw SireIO::parse_error(QObject::tr(
                                      "Cannot read the AmberRst file '%1' as it is missing the required 'Conventions' "
                                      "attribute.")
                                      .arg(this->filename()),
                                  CODELOC);
    }

    double convention_version;

    try
    {
        bool version_ok;

        convention_version = this->_lkr_getStringAttribute("ConventionVersion").toDouble(&version_ok);

        if (not version_ok)
        {
            // have to assume it is version 1...
            convention_version = 1;
        }
    }
    catch (...)
    {
        // have to assume it is version 1...
        convention_version = 1;
    }

    if (not conventions.contains("AMBER"))
    {
        // this is not an Amber RST file
        this->_lkr_close();
        throw SireIO::parse_error(QObject::tr(
                                      "Cannot read the AmberRst file '%1' as it doesn't match the required format. "
                                      "It is missing the 'AMBER' flag in the 'Conventions' attribute.")
                                      .arg(this->filename()),
                                  CODELOC);
    }

    is_restart = conventions.contains("AMBERRESTART");

    if (convention_version != 1)
    {
        this->_lkr_close();
        throw SireIO::parse_error(QObject::tr(
                                      "Cannot read the AmberRst file '%1' as it is written in version %1 of the format, "
                                      "and we only support version 1.0")
                                      .arg(convention_version),
                                  CODELOC);
    }

    try
    {
        application = this->_lkr_getStringAttribute("application");
    }
    catch (...)
    {
    }

    try
    {
        program = this->_lkr_getStringAttribute("program");
    }
    catch (...)
    {
    }

    try
    {
        program_version = this->_lkr_getStringAttribute("programVersion");
    }
    catch (...)
    {
    }

    try
    {
        title = this->_lkr_getStringAttribute("title");
    }
    catch (...)
    {
    }

    const auto dims = this->_lkr_getDimensions();

    nframes = dims.value("frame", 1);
    natoms = dims.value("atom", 0);

    int spatial = dims.value("spatial", 3);

    if (spatial != 3)
    {
        this->_lkr_close();
        throw SireIO::parse_error(QObject::tr(
                                      "We can only read AmberRst files with 3 spatial dimensions. The number "
                                      "of dimensions in '%1' is %2.")
                                      .arg(this->filename())
                                      .arg(spatial),
                                  CODELOC);
    }

    return true;
}

//////////
////////// Implementation of AmberRst
//////////

static const RegisterMetaType<AmberRst> r_rst;
const RegisterParser<AmberRst> register_rst;

QDataStream &operator<<(QDataStream &ds, const AmberRst &rst)
{
    writeHeader(ds, r_rst, 3);

    SharedDataStream sds(ds);

    sds << rst.current_frame << rst.parse_warnings
        << rst.nframes << rst.frame_idx
        << static_cast<const MoleculeParser &>(rst);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, AmberRst &rst)
{
    VersionID v = readHeader(ds, r_rst);

    if (v == 3)
    {
        SharedDataStream sds(ds);

        sds >> rst.current_frame >> rst.parse_warnings >> rst.nframes >> rst.frame_idx >> static_cast<MoleculeParser &>(rst);
    }
    else
        throw version_error(v, "3", r_rst, CODELOC);

    return ds;
}

/** Constructor */
AmberRst::AmberRst()
    : ConcreteProperty<AmberRst, MoleculeParser>(),
      nframes(0), frame_idx(0)
{
}

/** Return the format name that is used to identify this file format within Sire */
QString AmberRst::formatName() const
{
    return "RST";
}

/** Return the suffixes that AmberRst files will typically have */
QStringList AmberRst::formatSuffix() const
{
    static const QStringList suffixes = {"rst", "crd", "trj", "traj"};
    return suffixes;
}

/** Return a description of the file format */
QString AmberRst::formatDescription() const
{
    return QObject::tr("Amber coordinate/velocity binary (netcdf) restart/trajectory files "
                       "supported since Amber 9, now default since Amber 16.");
}

/** This is not a text file */
bool AmberRst::isTextFile() const
{
    return false;
}

/** Open the file and read in all the metadata */
void AmberRst::parse()
{
    // the netcdf conventions for Amber restart/trajectory files are
    // given here - http://ambermd.org/netcdf/nctraj.xhtml

    // NOTE THAT YOU CANNOT READ A NETCDF FILE IN PARALLEL
    // THE NETCDF LIBRARY IS NOT THREAD SAFE

    f.reset(new AmberRstFile(this->filename()));

    try
    {
        if (not f->open(QIODevice::ReadOnly, true, false))
        {
            throw SireIO::parse_error(QObject::tr(
                                          "Failed to open AmberRst file %1")
                                          .arg(this->filename()),
                                      CODELOC);
        }

        nframes = f->nFrames();

        if (nframes == 0)
        {
            this->setScore(0);
            f.reset();
            return;
        }

        current_frame = f->readFrame(0, this->usesParallel());
        frame_idx = 0;

        this->setScore(f->nFrames() * current_frame.nAtoms());
    }
    catch (...)
    {
        this->setScore(0);
        f.reset();
        throw;
    }
}

/** Construct by parsing the passed file */
AmberRst::AmberRst(const QString &filename, const PropertyMap &map)
    : ConcreteProperty<AmberRst, MoleculeParser>(map),
      nframes(0), frame_idx(0)
{
    // this gets the absolute file path
    this->setFilename(filename);
    this->parse();
}

/** Construct by parsing the data in the passed text lines */
AmberRst::AmberRst(const QStringList &lines, const PropertyMap &map)
    : ConcreteProperty<AmberRst, MoleculeParser>(lines, map)
{
    throw SireIO::parse_error(QObject::tr("You cannot create a binary AmberRst file from a set of text lines!"),
                              CODELOC);
}

/** Construct by extracting the necessary data from the passed System */
AmberRst::AmberRst(const System &system, const PropertyMap &map)
    : ConcreteProperty<AmberRst, MoleculeParser>(system, map),
      nframes(1), frame_idx(0)
{
    current_frame = MoleculeParser::createFrame(system, map);
}

/** Copy constructor */
AmberRst::AmberRst(const AmberRst &other)
    : ConcreteProperty<AmberRst, MoleculeParser>(other),
      current_frame(other.current_frame), parse_warnings(other.parse_warnings),
      nframes(other.nframes), frame_idx(other.frame_idx), f(other.f)
{
}

/** Destructor */
AmberRst::~AmberRst()
{
}

AmberRst &AmberRst::operator=(const AmberRst &other)
{
    if (this != &other)
    {
        current_frame = other.current_frame;
        parse_warnings = other.parse_warnings;
        nframes = other.nframes;
        frame_idx = other.frame_idx;
        f = other.f;

        MoleculeParser::operator=(other);
    }

    return *this;
}

bool AmberRst::operator==(const AmberRst &other) const
{
    return MoleculeParser::operator==(other);
}

bool AmberRst::operator!=(const AmberRst &other) const
{
    return MoleculeParser::operator!=(other);
}

const char *AmberRst::typeName()
{
    return QMetaType::typeName(qMetaTypeId<AmberRst>());
}

const char *AmberRst::what() const
{
    return AmberRst::typeName();
}

bool AmberRst::isFrame() const
{
    return true;
}

int AmberRst::nFrames() const
{
    return nframes;
}

int AmberRst::count() const
{
    return this->nFrames();
}

int AmberRst::size() const
{
    return this->nFrames();
}

Frame AmberRst::getFrame(int frame) const
{
    frame = SireID::Index(frame).map(this->nFrames());

    if (frame < 0)
        frame = 0;

    if (frame == frame_idx)
        return current_frame;

    if (f.get() == 0)
    {
        throw SireError::file_error(QObject::tr(
                                        "Somehow we don't have access to the underlying AmberRst file?"),
                                    CODELOC);
    }

    return f->readFrame(frame, this->usesParallel());
}

AmberRst AmberRst::operator[](int i) const
{
    i = SireID::Index(i).map(this->nFrames());

    AmberRst ret(*this);

    ret.current_frame = this->getFrame(i);
    ret.frame_idx = i;

    return ret;
}

QString AmberRst::toString() const
{
    if (this->nAtoms() == 0)
    {
        return QObject::tr("AmberRst::null");
    }
    else
    {
        return QObject::tr("AmberRst( nAtoms() = %1, nFrames() = %2 )")
            .arg(this->nAtoms())
            .arg(this->nFrames());
    }
}

/** Parse from the passed file */
AmberRst AmberRst::parse(const QString &filename)
{
    return AmberRst(filename);
}

/** Internal function used to add the data from this parser into the passed System */
void AmberRst::addToSystem(System &system, const PropertyMap &map) const
{
    MoleculeParser::copyFromFrame(current_frame, system, map);

    // update the System fileformat property to record that it includes
    // data from this file format
    QString fileformat = this->formatName();

    PropertyName fileformat_property = map["fileformat"];

    try
    {
        QString last_format = system.property(fileformat_property).asA<StringProperty>().value();
        fileformat = QString("%1,%2").arg(last_format, fileformat);
    }
    catch (...)
    {
    }

    if (fileformat_property.hasSource())
    {
        system.setProperty(fileformat_property.source(), StringProperty(fileformat));
    }
    else
    {
        system.setProperty("fileformat", StringProperty(fileformat));
    }
}

/** Return the number of atoms whose coordinates are contained in this restart file */
int AmberRst::nAtoms() const
{
    return current_frame.nAtoms();
}

/** Return this parser constructed from the passed filename */
MoleculeParserPtr AmberRst::construct(const QString &filename, const PropertyMap &map) const
{
    // don't construct from a pointer as it could leak
    return MoleculeParserPtr(AmberRst(filename, map));
}

/** Return this parser constructed from the passed set of lines */
MoleculeParserPtr AmberRst::construct(const QStringList &lines, const PropertyMap &map) const
{
    // don't construct from a pointer as it could leak
    return MoleculeParserPtr(AmberRst(lines, map));
}

/** Return this parser constructed from the passed SireSystem::System */
MoleculeParserPtr AmberRst::construct(const SireSystem::System &system, const PropertyMap &map) const
{
    // don't construct from a pointer as it could leak
    return MoleculeParserPtr(AmberRst(system, map));
}

/** Write this binary file 'filename' */
QStringList AmberRst::writeToFile(const QString &filename) const
{
    if (this->nFrames() == 0 or this->nAtoms() == 0)
        return QStringList();

    auto gil = SireBase::release_gil();

    createDirectoryForFile(filename);

    AmberRstFile outfile(filename);

    if (not outfile.open(QIODevice::WriteOnly))
        throw SireError::file_error(QObject::tr(
                                        "Could not open %1 to write the AmberRst file.")
                                        .arg(filename),
                                    CODELOC);

    if (this->writingTrajectory())
    {
        const auto frames = this->framesToWrite();

        ProgressBar bar("Save AmberRst", frames.count());
        bar.setSpeedUnit("frames / s");

        bar = bar.enter();

        for (int i = 0; i < frames.count(); ++i)
        {
            const auto frame = this->createFrame(frames[i]);
            outfile.writeFrame(frame, usesParallel());
            bar.setProgress(i + 1);
        }

        bar.success();
    }
    else
    {
        outfile.writeFrame(current_frame, usesParallel());
    }

    outfile.close();

    return QStringList(filename);
}
