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
  *  at https://sire.openbiosim.org
  *
\*********************************************/

#include "movingharmonicrestraints.h"

#include "SireUnits/units.h"

#include "SireID/index.h"

#include "SireError/errors.h"

#include "SireStream/datastream.h"
#include "SireStream/shareddatastream.h"

#include <QDebug>

using namespace SireMM;
using namespace SireMaths;
using namespace SireBase;
using namespace SireUnits;
using namespace SireUnits::Dimension;
using namespace SireStream;

///////
/////// Implementation of MovingHarmonicRestraint
///////

static const RegisterMetaType<MovingHarmonicRestraint> r_movingharmonicrest;

QDataStream &operator<<(QDataStream &ds, const MovingHarmonicRestraint &movingharmonicrest)
{
    writeHeader(ds, r_movingharmonicrest, 1);

    SharedDataStream sds(ds);

    sds << movingharmonicrest.atms0 << movingharmonicrest.atms1 << movingharmonicrest._k << movingharmonicrest._r0 << movingharmonicrest._r1
        << static_cast<const Property &>(movingharmonicrest);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, MovingHarmonicRestraint &movingharmonicrest)
{
    VersionID v = readHeader(ds, r_movingharmonicrest);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> movingharmonicrest.atms0 >> movingharmonicrest.atms1 >> movingharmonicrest._k >> movingharmonicrest._r0 >> movingharmonicrest._r1 >> static_cast<Property &>(movingharmonicrest);
    }
    else
        throw version_error(v, "1", r_movingharmonicrest, CODELOC);

    return ds;
}

/** Null constructor */
MovingHarmonicRestraint::MovingHarmonicRestraint()
    : ConcreteProperty<MovingHarmonicRestraint, Property>(),
      _k(0), _r0(0), _r1(0)
{
}

/** Construct to restrain the atom at index 'atom' to the specified position
 *  using the specified force constant and flat-bottom well-width
 */
MovingHarmonicRestraint::MovingHarmonicRestraint(qint64 atom0, qint64 atom1,
                             const SireUnits::Dimension::HarmonicBondConstant &k,
                             const SireUnits::Dimension::Length &r0,
                             const SireUnits::Dimension::Length &r1)
    : ConcreteProperty<MovingHarmonicRestraint, Property>(),
      _k(k), _r0(r0), _r1(r1)
{
    if (atom0 == atom1)
        throw SireError::invalid_arg(QObject::tr(
                                         "You cannot create a MovingHarmonic restraint between identical atoms! %1-%2")
                                         .arg(atom0)
                                         .arg(atom1),
                                     CODELOC);

    atms0 = QVector<qint64>(1, atom0);
    atms0.squeeze();

    atms1 = QVector<qint64>(1, atom1);
    atms1.squeeze();
}

/** Construct to restrain the centroid of the atoms whose indicies are
 *  in 'atoms' to the specified position using the specified force constant
 *  and flat-bottom well width
 */
MovingHarmonicRestraint::MovingHarmonicRestraint(const QList<qint64> &atoms0,
                             const QList<qint64> &atoms1,
                             const SireUnits::Dimension::HarmonicBondConstant &k,
                             const SireUnits::Dimension::Length &r0,
                             const SireUnits::Dimension::Length &r1)
    : ConcreteProperty<MovingHarmonicRestraint, Property>(),
      _k(k), _r0(r0), _r1(r1)
{
    if (atoms0.isEmpty() or atoms1.isEmpty())
        return;

    // remove duplicates
    atms0.reserve(atoms0.count());

    auto sorted = atoms0;
    std::sort(sorted.begin(), sorted.end());

    atms0.append(sorted.at(0));

    for (const auto &atom : sorted)
    {
        if (atom != atms0.last())
            atms0.append(atom);
    }

    atms0.squeeze();

    // now the same for atoms1
    atms1.reserve(atoms1.count());

    sorted = atoms1;
    std::sort(sorted.begin(), sorted.end());

    if (atms0.indexOf(sorted.at(0)) != -1)
        throw SireError::invalid_arg(QObject::tr(
                                         "You cannot have an overlap in atoms between the two groups. Atom "
                                         "%1 appears in both!")
                                         .arg(sorted.at(0)),
                                     CODELOC);

    atms1.append(sorted.at(0));

    for (const auto &atom : sorted)
    {
        if (atom != atms1.last())
        {
            if (atms0.indexOf(atom) != -1)
                throw SireError::invalid_arg(QObject::tr(
                                                 "You cannot have an overlap in atoms between the two groups. Atom "
                                                 "%1 appears in both!")
                                                 .arg(sorted.at(0)),
                                             CODELOC);

            atms1.append(atom);
        }

        atms1.squeeze();
    }
}

/** Copy constructor */
MovingHarmonicRestraint::MovingHarmonicRestraint(const MovingHarmonicRestraint &other)
    : ConcreteProperty<MovingHarmonicRestraint, Property>(other),
      atms0(other.atms0), atms1(other.atms1), _k(other._k), _r0(other._r0), _r1(other._r1)
{
}

MovingHarmonicRestraint::~MovingHarmonicRestraint()
{
}

MovingHarmonicRestraint &MovingHarmonicRestraint::operator=(const MovingHarmonicRestraint &other)
{
    if (this != &other)
    {
        atms0 = other.atms0;
        atms1 = other.atms1;
        _k = other._k;
        _r0 = other._r0;
        _r1 = other._r1;
    }

    return *this;
}

bool MovingHarmonicRestraint::operator==(const MovingHarmonicRestraint &other) const
{
    return atms0 == other.atms0 and atms1 == other.atms1 and
           _k == other._k and _r0 == other._r0 and _r1 == other._r1;
}

bool MovingHarmonicRestraint::operator!=(const MovingHarmonicRestraint &other) const
{
    return not operator==(other);
}

MovingHarmonicRestraints MovingHarmonicRestraint::operator+(const MovingHarmonicRestraint &other) const
{
    return MovingHarmonicRestraints(*this) + other;
}

MovingHarmonicRestraints MovingHarmonicRestraint::operator+(const MovingHarmonicRestraints &other) const
{
    return MovingHarmonicRestraints(*this) + other;
}

const char *MovingHarmonicRestraint::typeName()
{
    return QMetaType::typeName(qMetaTypeId<MovingHarmonicRestraint>());
}

const char *MovingHarmonicRestraint::what() const
{
    return MovingHarmonicRestraint::typeName();
}

MovingHarmonicRestraint *MovingHarmonicRestraint::clone() const
{
    return new MovingHarmonicRestraint(*this);
}

bool MovingHarmonicRestraint::isNull() const
{
    return atms0.isEmpty() or atms1.isEmpty();
}

QString MovingHarmonicRestraint::toString() const
{
    if (this->isNull())
        return QObject::tr("MovingHarmonicRestraint::null");

    else if (this->isAtomRestraint())
    {
        return QString("MovingHarmonicRestraint( %1 <=> %2, k=%3 : r0=%4 : r1=%5 )")
            .arg(this->atom0())
            .arg(this->atom1())
            .arg(_k.toString())
            .arg(_r0.toString())
            .arg(_r1.toString());
    }
    else
    {
        QStringList a0, a1;

        for (const auto &atom : atms0)
        {
            a0.append(QString::number(atom));
        }

        for (const auto &atom : atms1)
        {
            a1.append(QString::number(atom));
        }

        return QString("MovingHarmonicRestraint( [%1] <=> [%2], k=%3 : r0=%4 : r1=%5 )")
            .arg(a0.join(", "))
            .arg(a1.join(", "))
            .arg(_k.toString())
            .arg(_r0.toString())
            .arg(_r1.toString());
    }
}

/** Return whether this is a single-atom restraint */
bool MovingHarmonicRestraint::isAtomRestraint() const
{
    return atms0.count() == 1 and atms1.count() == 1;
}

/** Return whether this restraint acts on the centroid of a group
 *  of atoms */
bool MovingHarmonicRestraint::isCentroidRestraint() const
{
    return atms0.count() > 1 or atms1.count() > 1;
}

/** Return the index of the atom if this is a single-atom restraint */
qint64 MovingHarmonicRestraint::atom0() const
{
    if (not this->isAtomRestraint())
        throw SireError::incompatible_error(QObject::tr(
                                                "You cannot get the atom when this isn't a single-atom restraint!"),
                                            CODELOC);

    return this->atms0.at(0);
}

/** Return the index of the atom if this is a single-atom restraint */
qint64 MovingHarmonicRestraint::atom1() const
{
    if (not this->isAtomRestraint())
        throw SireError::incompatible_error(QObject::tr(
                                                "You cannot get the atom when this isn't a single-atom restraint!"),
                                            CODELOC);

    return this->atms1.at(0);
}

/** Return the indexes of the atoms whose centroid is to be restrained */
QVector<qint64> MovingHarmonicRestraint::atoms0() const
{
    if (not this->isCentroidRestraint())
        throw SireError::incompatible_error(QObject::tr(
                                                "You cannot get the atoms when this isn't a centroid restraint!"),
                                            CODELOC);

    return this->atms0;
}

/** Return the indexes of the atoms whose centroid is to be restrained */
QVector<qint64> MovingHarmonicRestraint::atoms1() const
{
    if (not this->isCentroidRestraint())
        throw SireError::incompatible_error(QObject::tr(
                                                "You cannot get the atoms when this isn't a centroid restraint!"),
                                            CODELOC);

    return this->atms1;
}

/** Return the force constant for the restraint */
SireUnits::Dimension::HarmonicBondConstant MovingHarmonicRestraint::k() const
{
    return this->_k;
}

/** Return the width of the harmonic MovingHarmonic. */
SireUnits::Dimension::Length MovingHarmonicRestraint::r0() const
{
    return this->_r0;
}

/** Return the dissociation energy (depth) of MovingHarmonic bond */
// SireUnits::Dimension::MolarEnergy MovingHarmonicRestraint::r1() const
SireUnits::Dimension::Length MovingHarmonicRestraint::r1() const
{
    return this->_r1;
}

///////
/////// Implementation of MovingHarmonicRestraints
///////

static const RegisterMetaType<MovingHarmonicRestraints> r_movingharmonicrests;

QDataStream &operator<<(QDataStream &ds, const MovingHarmonicRestraints &movingharmonicrests)
{
    writeHeader(ds, r_movingharmonicrests, 1);

    SharedDataStream sds(ds);

    sds << movingharmonicrests.r
        << static_cast<const Restraints &>(movingharmonicrests);

    return ds;
}

QDataStream &operator>>(QDataStream &ds, MovingHarmonicRestraints &movingharmonicrests)
{
    VersionID v = readHeader(ds, r_movingharmonicrests);

    if (v == 1)
    {
        SharedDataStream sds(ds);

        sds >> movingharmonicrests.r >> static_cast<Restraints &>(movingharmonicrests);
    }
    else
        throw version_error(v, "1", r_movingharmonicrests, CODELOC);

    return ds;
}

/** Null constructor */
MovingHarmonicRestraints::MovingHarmonicRestraints()
    : ConcreteProperty<MovingHarmonicRestraints, Restraints>()
{
}

MovingHarmonicRestraints::MovingHarmonicRestraints(const QString &name)
    : ConcreteProperty<MovingHarmonicRestraints, Restraints>(name)
{
}

MovingHarmonicRestraints::MovingHarmonicRestraints(const MovingHarmonicRestraint &restraint)
    : ConcreteProperty<MovingHarmonicRestraints, Restraints>()
{
    if (not restraint.isNull())
        r.append(restraint);
}

MovingHarmonicRestraints::MovingHarmonicRestraints(const QList<MovingHarmonicRestraint> &restraints)
    : ConcreteProperty<MovingHarmonicRestraints, Restraints>()
{
    for (const auto &restraint : restraints)
    {
        if (not restraint.isNull())
            r.append(restraint);
    }
}

MovingHarmonicRestraints::MovingHarmonicRestraints(const QString &name,
                               const MovingHarmonicRestraint &restraint)
    : ConcreteProperty<MovingHarmonicRestraints, Restraints>(name)
{
    if (not restraint.isNull())
        r.append(restraint);
}

MovingHarmonicRestraints::MovingHarmonicRestraints(const QString &name,
                               const QList<MovingHarmonicRestraint> &restraints)
    : ConcreteProperty<MovingHarmonicRestraints, Restraints>(name)
{
    for (const auto &restraint : restraints)
    {
        if (not restraint.isNull())
            r.append(restraint);
    }
}

MovingHarmonicRestraints::MovingHarmonicRestraints(const MovingHarmonicRestraints &other)
    : ConcreteProperty<MovingHarmonicRestraints, Restraints>(other), r(other.r)
{
}

MovingHarmonicRestraints::~MovingHarmonicRestraints()
{
}

MovingHarmonicRestraints &MovingHarmonicRestraints::operator=(const MovingHarmonicRestraints &other)
{
    r = other.r;
    Restraints::operator=(other);
    return *this;
}

bool MovingHarmonicRestraints::operator==(const MovingHarmonicRestraints &other) const
{
    return r == other.r and Restraints::operator==(other);
}

bool MovingHarmonicRestraints::operator!=(const MovingHarmonicRestraints &other) const
{
    return not operator==(other);
}

const char *MovingHarmonicRestraints::typeName()
{
    return QMetaType::typeName(qMetaTypeId<MovingHarmonicRestraints>());
}

const char *MovingHarmonicRestraints::what() const
{
    return MovingHarmonicRestraints::typeName();
}

MovingHarmonicRestraints *MovingHarmonicRestraints::clone() const
{
    return new MovingHarmonicRestraints(*this);
}

QString MovingHarmonicRestraints::toString() const
{
    if (this->isEmpty())
        return QObject::tr("MovingHarmonicRestraints::null");

    QStringList parts;

    const auto n = this->count();

    if (n <= 10)
    {
        for (int i = 0; i < n; ++i)
        {
            parts.append(QObject::tr("%1: %2").arg(i).arg(this->r.at(i).toString()));
        }
    }
    else
    {
        for (int i = 0; i < 5; ++i)
        {
            parts.append(QObject::tr("%1: %2").arg(i).arg(this->r.at(i).toString()));
        }

        parts.append("...");

        for (int i = n - 5; i < n; ++i)
        {
            parts.append(QObject::tr("%1: %2").arg(i).arg(this->r.at(i).toString()));
        }
    }

    return QObject::tr("MovingHarmonicRestraints( name=%1, size=%2\n%3\n)").arg(this->name()).arg(n).arg(parts.join("\n"));
}

/** Return whether or not this is empty */
bool MovingHarmonicRestraints::isEmpty() const
{
    return this->r.isEmpty();
}

/** Return whether or not this is empty */
bool MovingHarmonicRestraints::isNull() const
{
    return this->isEmpty();
}

/** Return the number of restraints */
int MovingHarmonicRestraints::nRestraints() const
{
    return this->r.count();
}

/** Return the number of restraints */
int MovingHarmonicRestraints::count() const
{
    return this->nRestraints();
}

/** Return the number of restraints */
int MovingHarmonicRestraints::size() const
{
    return this->nRestraints();
}

/** Return the number of atom restraints */
int MovingHarmonicRestraints::nAtomRestraints() const
{
    int n = 0;

    for (const auto &restraint : this->r)
    {
        n += int(restraint.isAtomRestraint());
    }

    return n;
}

/** Return the number of centroid restraints */
int MovingHarmonicRestraints::nCentroidRestraints() const
{
    int n = 0;

    for (const auto &restraint : this->r)
    {
        n += int(restraint.isCentroidRestraint());
    }

    return n;
}

/** Return whether or not there are any atom restraints */
bool MovingHarmonicRestraints::hasAtomRestraints() const
{
    for (const auto &restraint : this->r)
    {
        if (restraint.isAtomRestraint())
            return true;
    }

    return false;
}

/** Return whether or not there are any centroid restraints */
bool MovingHarmonicRestraints::hasCentroidRestraints() const
{
    for (const auto &restraint : this->r)
    {
        if (restraint.isCentroidRestraint())
            return true;
    }

    return false;
}

/** Return the ith restraint */
const MovingHarmonicRestraint &MovingHarmonicRestraints::at(int i) const
{
    i = SireID::Index(i).map(this->r.count());

    return this->r.at(i);
}

/** Return the ith restraint */
const MovingHarmonicRestraint &MovingHarmonicRestraints::operator[](int i) const
{
    return this->at(i);
}

/** Return all of the restraints */
QList<MovingHarmonicRestraint> MovingHarmonicRestraints::restraints() const
{
    return this->r;
}

/** Return all of the atom restraints */
QList<MovingHarmonicRestraint> MovingHarmonicRestraints::atomRestraints() const
{
    if (this->hasCentroidRestraints())
    {
        QList<MovingHarmonicRestraint> ar;

        for (const auto &restraint : this->r)
        {
            if (restraint.isAtomRestraint())
                ar.append(restraint);
        }

        return ar;
    }
    else
        return this->restraints();
}

/** Return all of the centroid restraints */
QList<MovingHarmonicRestraint> MovingHarmonicRestraints::centroidRestraints() const
{
    if (this->hasAtomRestraints())
    {
        QList<MovingHarmonicRestraint> cr;

        for (const auto &restraint : this->r)
        {
            if (restraint.isCentroidRestraint())
                cr.append(restraint);
        }

        return cr;
    }
    else
        return this->restraints();
}

/** Add a restraint onto the list */
void MovingHarmonicRestraints::add(const MovingHarmonicRestraint &restraint)
{
    if (not restraint.isNull())
        this->r.append(restraint);
}

/** Add a restraint onto the list */
void MovingHarmonicRestraints::add(const MovingHarmonicRestraints &restraints)
{
    this->r += restraints.r;
}

/** Add a restraint onto the list */
MovingHarmonicRestraints &MovingHarmonicRestraints::operator+=(const MovingHarmonicRestraint &restraint)
{
    this->add(restraint);
    return *this;
}

/** Add a restraint onto the list */
MovingHarmonicRestraints MovingHarmonicRestraints::operator+(const MovingHarmonicRestraint &restraint) const
{
    MovingHarmonicRestraints ret(*this);
    ret += restraint;
    return *this;
}

/** Add restraints onto the list */
MovingHarmonicRestraints &MovingHarmonicRestraints::operator+=(const MovingHarmonicRestraints &restraints)
{
    this->add(restraints);
    return *this;
}

/** Add restraints onto the list */
MovingHarmonicRestraints MovingHarmonicRestraints::operator+(const MovingHarmonicRestraints &restraints) const
{
    MovingHarmonicRestraints ret(*this);
    ret += restraints;
    return *this;
}