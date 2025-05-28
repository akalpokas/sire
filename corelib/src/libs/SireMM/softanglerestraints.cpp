/********************************************\
  *
  *  Sire - Molecular Simulation Framework
  *
  *  Copyright (C) 2025  Christopher Woods
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
  *  You can contact the authors at https://sire.openbiosim.org
  *
\*********************************************/

#include "softanglerestraints.h"

#include "SireID/index.h"

#include "SireUnits/units.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include "SireCAS/errors.h"

#include <QDebug>

using namespace SireMM;
using namespace SireID;
using namespace SireBase;
using namespace SireMaths;
using namespace SireStream;
using namespace SireUnits;
using namespace SireUnits::Dimension;

////////////
//////////// Implementation of SoftAngleRestraint
////////////

static const RegisterMetaType<SoftAngleRestraint> r_softangrest;

/** Serialise to a binary datastream */

QDataStream &operator<<(QDataStream &ds, const SoftAngleRestraint &softangrest)
{
    writeHeader(ds, r_softangrest, 1);

    SharedDataStream sds(ds);

    sds << softangrest.atms << softangrest._theta0 << softangrest._ktheta << softangrest._pe;

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, SoftAngleRestraint &softangrest)
{
    VersionID v = readHeader(ds, r_softangrest);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> softangrest.atms >> softangrest._theta0 >> softangrest._ktheta >> softangrest._pe;
    }
    else
        throw version_error(v, "1", r_softangrest, CODELOC);

    return ds;
}

/** Null constructor */
SoftAngleRestraint::SoftAngleRestraint()
    : ConcreteProperty<SoftAngleRestraint, Property>(),
      _ktheta(0), _theta0(0), _pe(0)
{
}

/** Construct a restraint that acts on the angle within the
    three atoms 'atom0', 'atom1' and 'atom2' (theta == a(012)),
    restraining the angle within these atoms */
SoftAngleRestraint::SoftAngleRestraint(const QList<qint64> &atoms,
                               const SireUnits::Dimension::Angle &theta0,
                               const SireUnits::Dimension::HarmonicAngleConstant &ktheta,
                               const SireUnits::Dimension::MolarEnergy &pe)
    : ConcreteProperty<SoftAngleRestraint, Property>(),
      _theta0(theta0), _ktheta(ktheta), _pe(pe)
{

    // Make sure that we have 3 distinct atoms
    QSet<qint64> distinct;
    distinct.reserve(3);

    for (const auto &atom : atoms)
    {
        if (atom >= 0)
            distinct.insert(atom);
    }

    atms = atoms.toVector();
}

/* Copy constructor*/
SoftAngleRestraint::SoftAngleRestraint(const SoftAngleRestraint &other)
    : ConcreteProperty<SoftAngleRestraint, Property>(other),
      atms(other.atms), _theta0(other._theta0), _ktheta(other._ktheta), _pe(other._pe)

{
}

/* Destructor */
SoftAngleRestraint::~SoftAngleRestraint()
{
}

SoftAngleRestraint &SoftAngleRestraint::operator=(const SoftAngleRestraint &other)
{
    if (this != &other)
    {
        Property::operator=(other);
        atms = other.atms;
        _theta0 = other._theta0;
        _ktheta = other._ktheta;
        _pe = other._pe;
    }

    return *this;
}

bool SoftAngleRestraint::operator==(const SoftAngleRestraint &other) const
{
    return atms == other.atms and
           _theta0 == other._theta0 and
           _ktheta == other._ktheta and
           _pe == other._pe;
}

bool SoftAngleRestraint::operator!=(const SoftAngleRestraint &other) const
{
    return not operator==(other);
}

SoftAngleRestraints SoftAngleRestraint::operator+(const SoftAngleRestraint &other) const
{
    return SoftAngleRestraints(*this) + other;
}

SoftAngleRestraints SoftAngleRestraint::operator+(const SoftAngleRestraints &other) const
{
    return SoftAngleRestraints(*this) + other;
}

const char *SoftAngleRestraint::typeName()
{
    return QMetaType::typeName(qMetaTypeId<SoftAngleRestraint>());
}

const char *SoftAngleRestraint::what() const
{
    return SoftAngleRestraint::typeName();
}

SoftAngleRestraint *SoftAngleRestraint::clone() const
{
    return new SoftAngleRestraint(*this);
}

bool SoftAngleRestraint::isNull() const
{
    return atms.isEmpty();
}

QString SoftAngleRestraint::toString() const
{
    if (this->isNull())
        return QObject::tr("SoftAngleRestraint::null");
    else
    {
        QStringList a;

        for (const auto &atom : atms)
        {
            a.append(QString::number(atom));
        }
        return QString("SoftAngleRestraint( [%1], theta0=%2, ktheta=%3, pe=%4 )")
            .arg(a.join(", "))
            .arg(_theta0.toString())
            .arg(_ktheta.toString())
            .arg(_pe.toString());
    }
}

/** Return the force constant for the restraint */
SireUnits::Dimension::HarmonicAngleConstant SoftAngleRestraint::ktheta() const
{
    return this->_ktheta;
}

/** Return the equilibrium angle for the restraint */
SireUnits::Dimension::Angle SoftAngleRestraint::theta0() const
{
    return this->_theta0;
}

/** Return the plateau energy of the restraint */
SireUnits::Dimension::MolarEnergy SoftAngleRestraint::pe() const
{
    return this->_pe;
}

/** Return the atoms involved in the restraint */
QVector<qint64> SoftAngleRestraint::atoms() const
{
    return this->atms;
}

///////
/////// Implementation of SoftAngleRestraints
///////

/** Serialise to a binary datastream */

static const RegisterMetaType<SoftAngleRestraints> r_softangrests;

QDataStream &operator<<(QDataStream &ds, const SoftAngleRestraints &softangrests)
{
    writeHeader(ds, r_softangrests, 1);

    SharedDataStream sds(ds);

    sds << softangrests.r
        << static_cast<const Restraints &>(softangrests);

    return ds;
}

/** Extract from a binary datastream */
QDataStream &operator>>(QDataStream &ds, SoftAngleRestraints &softangrests)
{
    VersionID v = readHeader(ds, r_softangrests);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> softangrests.r >>
            static_cast<Restraints &>(softangrests);
    }
    else
        throw version_error(v, "1", r_softangrests, CODELOC);

    return ds;
}

/** Null constructor */
SoftAngleRestraints::SoftAngleRestraints()
    : ConcreteProperty<SoftAngleRestraints, Restraints>()
{
}

SoftAngleRestraints::SoftAngleRestraints(const QString &name)
    : ConcreteProperty<SoftAngleRestraints, Restraints>(name)
{
}

SoftAngleRestraints::SoftAngleRestraints(const SoftAngleRestraint &restraint)
    : ConcreteProperty<SoftAngleRestraints, Restraints>()
{
    if (not restraint.isNull())
        r.append(restraint);
}

SoftAngleRestraints::SoftAngleRestraints(const QList<SoftAngleRestraint> &restraints)
    : ConcreteProperty<SoftAngleRestraints, Restraints>()
{
    for (const auto &restraint : restraints)
    {
        if (not restraint.isNull())
            r.append(restraint);
    }
}

SoftAngleRestraints::SoftAngleRestraints(const QString &name,
                                 const SoftAngleRestraint &restraint)
    : ConcreteProperty<SoftAngleRestraints, Restraints>(name)
{
    if (not restraint.isNull())
        r.append(restraint);
}

SoftAngleRestraints::SoftAngleRestraints(const QString &name,
                                 const QList<SoftAngleRestraint> &restraints)
    : ConcreteProperty<SoftAngleRestraints, Restraints>(name)
{
    for (const auto &restraint : restraints)
    {
        if (not restraint.isNull())
            r.append(restraint);
    }
}

SoftAngleRestraints::SoftAngleRestraints(const SoftAngleRestraints &other)
    : ConcreteProperty<SoftAngleRestraints, Restraints>(other), r(other.r)
{
}

/* Desctructor */
SoftAngleRestraints::~SoftAngleRestraints()
{
}

SoftAngleRestraints &SoftAngleRestraints::operator=(const SoftAngleRestraints &other)
{
    if (this != &other)
    {
        Restraints::operator=(other);
        r = other.r;
    }

    return *this;
}

bool SoftAngleRestraints::operator==(const SoftAngleRestraints &other) const
{
    return r == other.r and Restraints::operator==(other);
}

bool SoftAngleRestraints::operator!=(const SoftAngleRestraints &other) const
{
    return not operator==(other);
}

const char *SoftAngleRestraints::typeName()
{
    return QMetaType::typeName(qMetaTypeId<SoftAngleRestraints>());
}

const char *SoftAngleRestraints::what() const
{
    return SoftAngleRestraints::typeName();
}

SoftAngleRestraints *SoftAngleRestraints::clone() const
{
    return new SoftAngleRestraints(*this);
}

QString SoftAngleRestraints::toString() const
{
    if (this->isEmpty())
        return QObject::tr("SoftAngleRestraints::null");

    QStringList parts;

    const auto n = this->count();

    if (n <= 10)
    {
        for (int i = 0; i < n; i++)
        {
            parts.append(QObject::tr("%1: %2").arg(i).arg(this->r.at(i).toString()));
        }
    }
    else
    {
        for (int i = 0; i < 5; i++)
        {
            parts.append(QObject::tr("%1: %2").arg(i).arg(this->r.at(i).toString()));
        }

        parts.append("...");

        for (int i = n - 5; i < n; i++)
        {
            parts.append(QObject::tr("%1: %2").arg(i).arg(this->r.at(i).toString()));
        }
    }

    return QObject::tr("SoftAngleRestraints( name=%1, size=%2\n%3\n )")
        .arg(this->name())
        .arg(n)
        .arg(parts.join("\n"));
}

/** Return whether or not this is empty */
bool SoftAngleRestraints::isEmpty() const
{
    return this->r.isEmpty();
}

/** Return whether or not this is empty */
bool SoftAngleRestraints::isNull() const
{
    return this->isEmpty();
}

/** Return the number of restraints */
int SoftAngleRestraints::nRestraints() const
{
    return this->r.count();
}

/** Return the number of restraints */
int SoftAngleRestraints::count() const
{
    return this->nRestraints();
}

/** Return the number of restraints */
int SoftAngleRestraints::size() const
{
    return this->nRestraints();
}

/** Return the ith restraint */
const SoftAngleRestraint &SoftAngleRestraints::at(int i) const
{
    i = SireID::Index(i).map(this->r.count());

    return this->r.at(i);
}

/** Return the ith restraint */
const SoftAngleRestraint &SoftAngleRestraints::operator[](int i) const
{
    return this->at(i);
}

/** Return all of the restraints */
QList<SoftAngleRestraint> SoftAngleRestraints::restraints() const
{
    return this->r;
}

/** Add a restraints onto the list */
void SoftAngleRestraints::add(const SoftAngleRestraint &restraint)
{
    if (not restraint.isNull())
        r.append(restraint);
}

/** Add a restraint onto the list */
void SoftAngleRestraints::add(const SoftAngleRestraints &restraints)
{
    this->r += restraints.r;
}

/** Add a restraint onto the list */
SoftAngleRestraints &SoftAngleRestraints::operator+=(const SoftAngleRestraint &restraint)
{
    this->add(restraint);
    return *this;
}

/** Add a restraint onto the list */
SoftAngleRestraints SoftAngleRestraints::operator+(const SoftAngleRestraint &restraint) const
{
    SoftAngleRestraints ret(*this);
    ret += restraint;
    return *this;
}

/** Add restraints onto the list */
SoftAngleRestraints &SoftAngleRestraints::operator+=(const SoftAngleRestraints &restraints)
{
    this->add(restraints);
    return *this;
}

/** Add restraints onto the list */
SoftAngleRestraints SoftAngleRestraints::operator+(const SoftAngleRestraints &restraints) const
{
    SoftAngleRestraints ret(*this);
    ret += restraints;
    return *this;
}
