/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2022  Christopher Woods
  *
  *  This program is free software; you can redistribute it and/or modify
  *  it under the terms of the GNU General Public License as published by
  *  the Free Software Foundation; either version 2 of the License, or
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
  *  You can contact the authors via the developer's mailing list
  *  at http://siremol.org
  *
\*********************************************/


#include "SireIO/sdf.h"

#include "SireSystem/system.h"

#include "SireBase/parallel.h"
#include "SireBase/stringproperty.h"

#include "SireError/errors.h"
#include "SireIO/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireMol/atomcharges.h"
#include "SireMol/atomcoords.h"
#include "SireMol/atomelements.h"
#include "SireMol/errors.h"
#include "SireMol/molecule.h"
#include "SireMol/moleditor.h"

#include "SireUnits/units.h"

#include <QFile>
#include <QtMath>

using namespace SireBase;
using namespace SireIO;
using namespace SireMol;
using namespace SireStream;
using namespace SireSystem;
using namespace SireUnits;

namespace SireIO { namespace detail {

class SDFAtom
{
public:
    SDFAtom()
    {}

    ~SDFAtom()
    {}

    QString toString() const
    {
        return QString("%1  %2 %3 %4  %5 %6  %7")
                .arg(name).arg(x).arg(y).arg(z)
                .arg(mass_difference).arg(chg_difference)
                .arg(fields.join(":"));
    }

    QString name;
    double x;
    double y;
    double z;
    qint32 mass_difference;
    qint32 chg_difference;
    QStringList fields;
};

class SDFBond
{
public:
    SDFBond()
    {}

    ~SDFBond()
    {}

    QString toString() const
    {
        return QString("%1-%2  %3  %4  %5")
                    .arg(atom0).arg(atom1)
                    .arg(typ).arg(stereoscopy)
                    .arg(fields.join(":"));
    }

    qint32 atom0;
    qint32 atom1;
    qint32 typ;
    qint32 stereoscopy;
    QStringList fields;
};

class SDFMolecule
{
public:
    SDFMolecule()
    {}

    ~SDFMolecule()
    {}

    QString toString() const
    {
        QStringList lines;

        lines.append(QString("name = %1").arg(name));
        lines.append(QString("software = %1").arg(software));
        lines.append(QString("comment = %1").arg(comment));

        lines.append(QString("counts: %1").arg(counts.join(" ")));

        lines.append(QString("nAtoms == %1").arg(atoms.count()));

        for (auto atom : atoms)
        {
            lines.append(atom.toString());
        }

        lines.append(QString("nBonds == %1").arg(bonds.count()));

        for (auto bond : bonds)
        {
            lines.append(bond.toString());
        }

        for (auto key : properties.keys())
        {
            lines.append(QString("property %1").arg(key));

            for (auto s : properties[key])
            {
                lines.append(s);
            }
        }

        for (auto key : data.keys())
        {
            lines.append(QString("data %1").arg(key));

            for (auto s : data[key])
            {
                lines.append(s);
            }
        }

        return lines.join("\n");
    }

    QString name;
    QString software;
    QString comment;
    QStringList counts;

    QVector<SDFAtom> atoms;
    QVector<SDFBond> bonds;

    QHash<QString, QStringList> properties;
    QHash<QString, QStringList> data;
};

}} // end of namespace SireIO::detail

QDataStream &operator<<(QDataStream &ds, const SireIO::detail::SDFAtom &atom)
{
    ds << atom.name << atom.x << atom.y << atom.z
       << atom.chg_difference << atom.mass_difference
       << atom.fields;

    return ds;
}

QDataStream &operator>>(QDataStream &ds, SireIO::detail::SDFAtom &atom)
{
    ds >> atom.name >> atom.x >> atom.y >> atom.z
       >> atom.chg_difference >> atom.mass_difference
       >> atom.fields;

    return ds;
}

QDataStream &operator<<(QDataStream &ds, const SireIO::detail::SDFBond &bond)
{
    ds << bond.atom0 << bond.atom1 << bond.typ
       << bond.stereoscopy << bond.fields;

    return ds;
}

QDataStream &operator>>(QDataStream &ds, SireIO::detail::SDFBond &bond)
{
    ds >> bond.atom0 >> bond.atom1 >> bond.typ
       >> bond.stereoscopy >> bond.fields;

    return ds;
}

QDataStream &operator<<(QDataStream &ds, const SireIO::detail::SDFMolecule &mol)
{
    ds << mol.atoms << mol.bonds << mol.properties << mol.data;
    return ds;
}

QDataStream &operator>>(QDataStream &ds, SireIO::detail::SDFMolecule &mol)
{
    ds >> mol.atoms >> mol.bonds >> mol.properties >> mol.data;
    return ds;
}

using namespace SireIO::detail;

const RegisterParser<SDF> register_sdf;
static const RegisterMetaType<SDF> r_sdf;

QDataStream &operator<<(QDataStream &ds, const SDF &sdf)
{
    writeHeader(ds, r_sdf, 1);

    SharedDataStream sds(ds);

    sds << sdf.molecules << sdf.parse_warnings
        << static_cast<const MoleculeParser&>(sdf);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, SDF &sdf)
{
    VersionID v = readHeader(ds, r_sdf);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> sdf.molecules >> sdf.parse_warnings
            >> static_cast<MoleculeParser&>(sdf);
    }
    else
        throw version_error(v, "1", r_sdf, CODELOC);

    return ds;
}

/** Constructor */
SDF::SDF() : ConcreteProperty<SDF,MoleculeParser>()
{}

/** Construct to read in the data from the file called 'filename'. The
    passed property map can be used to pass extra parameters to control
    the parsing */
SDF::SDF(const QString &filename, const PropertyMap &map) :
    ConcreteProperty<SDF,MoleculeParser>(filename,map)
{
    //the file has been read into memory and is available via
    //the MoleculeParser::lines() function.

    //a parameter has also been read in MoleculeParser to say whether
    //we are allowed to use multiple cores to parse the file, e.g.
    //MoleculeParser::usesParallel() will be true

    //parse the data in the parse function
    this->parseLines(map);

    //now make sure that everything is correct with this object
    this->assertSane();
}

/** Construct to read in the data from the passed text lines. The
    passed property map can be used to pass extra parameters to control
    the parsing */
SDF::SDF(const QStringList &lines, const PropertyMap &map) :
    ConcreteProperty<SDF,MoleculeParser>(lines,map)
{
    //the file has been read into memory and is available via
    //the MoleculeParser::lines() function.

    //a parameter has also been read in MoleculeParser to say whether
    //we are allowed to use multiple cores to parse the file, e.g.
    //MoleculeParser::usesParallel() will be true

    //parse the data in the parse function
    this->parseLines(map);

    //now make sure that everything is correct with this object
    this->assertSane();
}

/** Construct this parser by extracting all necessary information from the
    passed SireSystem::System, looking for the properties that are specified
    in the passed property map */
SDF::SDF(const SireSystem::System &system, const PropertyMap &map) :
    ConcreteProperty<SDF,MoleculeParser>(map)
{

    // write the file into lines, then parse them back for self-consistency
    QStringList lines;

    // Reparse the lines as a self-consistency check.
    SDF parsed(lines, map);

    this->operator=(parsed);
}

/** Copy constructor */
SDF::SDF(const SDF &other) :
    ConcreteProperty<SDF,MoleculeParser>(other),
    molecules(other.molecules),
    parse_warnings(other.parse_warnings)
{}

/** Destructor */
SDF::~SDF()
{}

/** Copy assignment operator */
SDF& SDF::operator=(const SDF &other)
{
    if (this != &other)
    {
        this->molecules = other.molecules;
        this->parse_warnings = other.parse_warnings;

        MoleculeParser::operator=(other);
    }

    return *this;
}

/** Comparison operator */
bool SDF::operator==(const SDF &other) const
{
    return MoleculeParser::operator==(other);
}

/** Comparison operator */
bool SDF::operator!=(const SDF &other) const
{
    return not operator==(other);
}

/** Return the C++ name for this class */
const char* SDF::typeName()
{
    return QMetaType::typeName( qMetaTypeId<SDF>() );
}

/** Return the C++ name for this class */
const char* SDF::what() const
{
    return SDF::typeName();
}

/** Return the parser that has been constructed by reading in the passed
    file using the passed properties */
MoleculeParserPtr SDF::construct(const QString &filename,
                                 const PropertyMap &map) const
{
    return SDF(filename,map);
}

/** Return the parser that has been constructed by reading in the passed
    text lines using the passed properties */
MoleculeParserPtr SDF::construct(const QStringList &lines,
                                 const PropertyMap &map) const
{
    return SDF(lines,map);
}

/** Return the parser that has been constructed by extract all necessary
    data from the passed SireSystem::System using the specified properties */
MoleculeParserPtr SDF::construct(const SireSystem::System &system,
                                 const PropertyMap &map) const
{
    return SDF(system,map);
}

/** Return a string representation of this parser */
QString SDF::toString() const
{
    if (lines().isEmpty())
        return QObject::tr("SDF::null");
    else
    {
        return QObject::tr("SDF( nMolecules() == %1 )")
                        .arg(this->nMolecules());
    }
}

/** Return the format name that is used to identify this file format within Sire */
QString SDF::formatName() const
{
    return "SDF";
}

/** Return any warnings raised when parsing this file */
QStringList SDF::parseWarnings() const
{
    return this->parse_warnings;
}

/** Return the number of molecules loaded in this file */
int SDF::nMolecules() const
{
    return this->molecules.count();
}

/** Return a description of the file format */
QString SDF::formatDescription() const
{
    return QObject::tr("Structure Data File (SDF) format files.");
}

/** Return the suffixes that these files are normally associated with */
QStringList SDF::formatSuffix() const
{
    static const QStringList suffixes = { "SDF" };
    return suffixes;
}

/** Return whether or not this is a lead parser. The lead parser is responsible
    for starting the process of turning the parsed file into the System. There
    must be one and one-only lead parser in a set of parsers creating a System */
bool SDF::isLead() const
{
    return true;
}

/** Function that is called to assert that this object is sane. This
    should raise an exception if the parser is in an invalid state */
void SDF::assertSane() const
{
    if (this->nMolecules() == 0)
    {
        QStringList errors = this->parse_warnings;

        if (not errors.isEmpty())
        {
            throw SireIO::parse_error(QObject::tr(
                "There were errors reading the SDF format "
                "file:\n%1").arg(errors.join("\n\n")), CODELOC);
        }
    }
}

/** Internal function that is used to actually parse the data contained
    in the lines of the file */
void SDF::parseLines(const PropertyMap &map)
{
    /* File format is decribed here:
        https://www.herongyang.com/Molecule/SDF-Format-Specification.html
        http://www.nonlinear.com/progenesis/sdf-studio/v0.9/faq/sdf-file-format-guidance.aspx
     */

    // The first three lines are the header block
    const auto &l = this->lines();

    if (l.count() < 4)
    {
        // this is not a valid SDF file
        this->parse_warnings.append(
            QObject::tr("There are no enough lines for this "
                        "to be a valid SDF-formatted file.")
        );
        this->setScore(0);
        return;
    }

    // first line is the molecule name
    const auto molname = l[0].simplified();
    // then the software used to generate the SDF file
    const auto software = l[1].simplified();
    // then a user-supplied comment (if any)
    const auto comment = l[2].simplified();

    // now the counts line - 12 fixed-width fields. The first eleven are 3
    // characters long. The last is 6 characters long
    const auto counts_line = l[3];

    if (counts_line.size() < 39)
    {
        this->parse_warnings.append(
            QObject::tr("The counts line in this SDF file does not "
                        "have enough characters! '%1'. It should be "
                        "at least 39 characters wide.").arg(counts_line)
        );
        this->setScore(0);
        return;
    }

    QStringList counts;

    for (int i=0; i<11; ++i)
    {
        counts.append(counts_line.mid(i*3, 3));
    }

    // the last counts line item can also be a string!
    counts.append(counts_line.mid(33, 6));

    bool ok;

    // Atom counter.
    const int natoms = counts[0].toInt(&ok);

    if (not ok)
    {
        this->setScore(0);
        this->parse_warnings.append(
            QObject::tr("Cannot interpret the number of atoms from the "
                        "counts line: %1").arg(counts_line));
    }

    // Bonds counter
    const int nbonds = counts[1].toInt(&ok);

    if (not ok)
    {
        this->setScore(0);
        this->parse_warnings.append(
            QObject::tr("Cannot interpret the number of bonds from the "
                        "counts line: %1").arg(counts_line));
    }

    if (natoms == 0)
    {
        // nothing to read?
        this->setScore(0);
        this->parse_warnings.append(
            QObject::tr("The number of atoms to read is set to zero?")
        );
        return;
    }

    // next, read in the atoms
    if (l.count() < 4 + natoms + nbonds)
    {
        this->parse_warnings.append(
            QObject::tr("There aren't enough lines in this file to "
                        "contain all of the atoms and bonds. File is "
                        "corrupted?")
        );
        this->setScore(0);
        return;
    }

    QVector<SDFAtom> atoms(natoms);

    for (int i=0; i<natoms; ++i)
    {
        QString line = l[i+4];

        //000000000011111111112222222222333333333344444444445555555555666666666
        //012345678901234567890123456789012345678901234567890123456789012345678
        //    0.5369    0.9749    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0
        while (line.size() < 69)
        {
            line.append(" ");
        }

        SDFAtom &atom = atoms[i];

        auto assert_ok = [&](bool is_ok, int line_num,
                             const QString &line, const QString &field)
                        {
                            if (not is_ok)
                            {
                                this->parse_warnings.append(
                                    QObject::tr("Atom line %1 has a problem "
                                      "with field '%2'. '%3'")
                                        .arg(line_num)
                                        .arg(field)
                                        .arg(line)
                                );

                                this->setScore(0);
                                return false;
                            }

                            return true;
                        };

        atom.x = line.midRef(0,10).toDouble(&ok);

        if (not assert_ok(ok, i+1, line, "x"))
            return;

        atom.y = line.midRef(10, 10).toDouble(&ok);

        if (not assert_ok(ok, i+1, line, "y"))
            return;

        atom.z = line.midRef(20, 10).toDouble(&ok);

        if (not assert_ok(ok, i+1, line, "z"))
            return;

        atom.name = line.mid(31,3);

        atom.mass_difference = line.midRef(34, 2).toInt(&ok);

        if (not assert_ok(ok, i+1, line, "mass difference"))
            return;

        if (atom.mass_difference < -3 or atom.mass_difference > 4)
        {
            this->parse_warnings.append(QObject::tr(
                "Only mass differences between -3 and 4 are supported. "
                "Cannot have a difference of %1 on line %2. '%3'")
                    .arg(atom.mass_difference)
                    .arg(i+1).arg(line));
            this->setScore(0);
            return;
        }

        atom.chg_difference = line.midRef(36, 3).toInt(&ok);

        if (not assert_ok(ok, i+1, line, "charge difference"))
            return;

        if (atom.chg_difference < 0 or atom.chg_difference > 7)
        {
            this->parse_warnings.append(QObject::tr(
                "Only charge differences between 0 and 7 are supported. "
                "Cannot have a difference of %1 on line %2. '%3'")
                    .arg(atom.chg_difference)
                    .arg(i+1).arg(line));
        }

        // ten more fields of 3 characters each. We won't convert these
        // to numbers - just leave as strings
        for (int j=0; j<10; ++j)
        {
            atom.fields.append(line.mid(38+(3*j), 3));
        }
    }

    QVector<SDFBond> bonds(nbonds);

    for (int i=0; i<nbonds; ++i)
    {
        QString line = l[4 + natoms + i];

        //000000000011111111112
        //012345678901234567890
        //  1  2  1  0  0  0  0
        while (line.size() < 20)
        {
            line.append(" ");
        }

        SDFBond &bond = bonds[i];

        auto assert_ok = [&](bool is_ok, int line_num,
                             const QString &line, const QString &field)
                        {
                            if (not is_ok)
                            {
                                this->parse_warnings.append(
                                    QObject::tr("Bond line %1 has a problem "
                                      "with field '%2'. '%3'")
                                        .arg(line_num)
                                        .arg(field)
                                        .arg(line)
                                );

                                this->setScore(0);
                                return false;
                            }

                            return true;
                        };

        bond.atom0 = line.midRef(0, 3).toInt(&ok);

        if (not assert_ok(ok, i+1, line, "atom0"))
            return;

        bond.atom1 = line.midRef(3, 3).toInt(&ok);

        if (not assert_ok(ok, i+1, line, "atom1"))
            return;

        bond.typ = line.midRef(6, 3).toInt(&ok);

        if (not assert_ok(ok, i+1, line, "bond type"))
            return;

        bond.stereoscopy = line.midRef(9, 3).toInt(&ok);

        if (not assert_ok(ok, i+1, line, "stereoscopy"))
            return;

        // now add on the three fields, which we will leave as strings
        for (int j=0; j<3; ++j)
        {
            bond.fields.append(line.mid(12+(3*j), 3));
        }
    }

    // now read in the properties...
    int linenum = 4 + natoms + nbonds;

    QHash<QString, QStringList> properties;

    while (linenum < l.count())
    {
        const auto line = l[linenum];

        if (line.startsWith("M  END"))
        {
            // end of the properties
            break;
        }
        else if (line.startsWith("M  "))
        {
            // this is a new property line
            if (line.size() >= 6)
            {
                QString key = line.mid(3, 3);
                properties[key].append(line);
            }
        }

        linenum += 1;
    }

    // now read in the data...
    QString key = "";
    QHash<QString, QStringList> data;

    while (linenum < l.count())
    {
        QString line = l[linenum];

        if (line.startsWith("$$$$"))
        {
            // end of the data (and molecule)
            break;
        }
        else if (line.startsWith("> "))
        {
            line = line.mid(2);

            int start_idx = line.indexOf("<");
            int end_idx = line.indexOf(">");

            if (start_idx >= 0 and end_idx >= 0)
            {
                key = line.mid(start_idx, end_idx-start_idx+1);
            }
            else
            {
                key = "";
            }
        }
        else if (key.size() > 0)
        {
            if (line.simplified().size() > 0)
            {
                data[key].append(line);
            }
        }

        linenum += 1;
    }

    SDFMolecule molecule;
    molecule.name = molname;
    molecule.software = software;
    molecule.comment = comment;
    molecule.counts = counts;
    molecule.atoms = atoms;
    molecule.bonds = bonds;
    molecule.properties = properties;
    molecule.data = data;

    auto debug_lines = molecule.toString().split("\n");

    for (auto l: debug_lines)
    {
        qDebug() << l;
    }

    this->molecules.append(molecule);

    this->setScore(1000);
}

/** Use the data contained in this parser to create a new System of molecules,
    assigning properties based on the mapping in 'map'. */
System SDF::startSystem(const PropertyMap &map) const
{
    return System();
}

/** Use the data contained in this parser to add information from the file to
    the molecules that exist already in the passed System. For example, this
    may be used to add coordinate data from this file to the molecules in
    the passed System that are missing coordinate data. */
void SDF::addToSystem(System &system, const PropertyMap &map) const
{
    /* Here we add information to an existing system that has been
       generated by a lead parser.
     */
    return;

  	// Update the System fileformat property to record that it includes
    // data from this file format.
    QString fileformat = this->formatName();

    PropertyName fileformat_property = map["fileformat"];

    try
    {
        QString last_format = system.property(fileformat_property)
                                    .asA<StringProperty>().value();

        fileformat = QString("%1,%2").arg(last_format,fileformat);
    }
    catch(...)
    {}

    if (fileformat_property.hasSource())
    {
        system.setProperty(fileformat_property.source(), StringProperty(fileformat));
    }
    else
    {
        system.setProperty("fileformat", StringProperty(fileformat));
    }
}
