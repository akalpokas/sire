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

#ifndef SIREMM_MOVINGHARMONICRESTRAINTS_H
#define SIREMM_MOVINGHARMONICRESTRAINTS_H

#include "restraints.h"

#include "SireUnits/dimensions.h"
#include "SireUnits/generalunit.h"

SIRE_BEGIN_HEADER

namespace SireMM
{
    class MovingHarmonicRestraint;
    class MovingHarmonicRestraints;
}

SIREMM_EXPORT QDataStream &operator<<(QDataStream &, const SireMM::MovingHarmonicRestraint &);
SIREMM_EXPORT QDataStream &operator>>(QDataStream &, SireMM::MovingHarmonicRestraint &);

SIREMM_EXPORT QDataStream &operator<<(QDataStream &, const SireMM::MovingHarmonicRestraints &);
SIREMM_EXPORT QDataStream &operator>>(QDataStream &, SireMM::MovingHarmonicRestraints &);

namespace SireMM
{

    /** This class represents a single Morse restraint between any two
     *  atoms in a system (or between the centroids of any two groups
     *  of atoms in a system)
     */
    class SIREMM_EXPORT MovingHarmonicRestraint
        : public SireBase::ConcreteProperty<MovingHarmonicRestraint, SireBase::Property>
    {
        friend QDataStream & ::operator<<(QDataStream &, const SireMM::MovingHarmonicRestraint &);
        friend QDataStream & ::operator>>(QDataStream &, SireMM::MovingHarmonicRestraint &);

    public:
        MovingHarmonicRestraint();
        MovingHarmonicRestraint(qint64 atom0, qint64 atom1,
                       const SireUnits::Dimension::HarmonicBondConstant &k,
                       const SireUnits::Dimension::Length &r0,
                       const SireUnits::Dimension::Length &r1);

        MovingHarmonicRestraint(const QList<qint64> &atoms0,
                       const QList<qint64> &atoms1,
                       const SireUnits::Dimension::HarmonicBondConstant &k,
                       const SireUnits::Dimension::Length &r0,
                       const SireUnits::Dimension::Length &r1);

        MovingHarmonicRestraint(const MovingHarmonicRestraint &other);

        ~MovingHarmonicRestraint();

        MovingHarmonicRestraint &operator=(const MovingHarmonicRestraint &other);

        bool operator==(const MovingHarmonicRestraint &other) const;
        bool operator!=(const MovingHarmonicRestraint &other) const;

        MovingHarmonicRestraints operator+(const MovingHarmonicRestraint &other) const;
        MovingHarmonicRestraints operator+(const MovingHarmonicRestraints &other) const;

        static const char *typeName();
        const char *what() const;

        MovingHarmonicRestraint *clone() const;

        QString toString() const;

        bool isNull() const;

        bool isAtomRestraint() const;
        bool isCentroidRestraint() const;

        qint64 atom0() const;
        qint64 atom1() const;

        QVector<qint64> atoms0() const;
        QVector<qint64> atoms1() const;

        SireUnits::Dimension::HarmonicBondConstant k() const;
        SireUnits::Dimension::Length r0() const;
        SireUnits::Dimension::Length r1() const;

    private:
        /** The first set of atoms involved in the restraint */
        QVector<qint64> atms0;

        /** The second set of atoms involved in the restraint */
        QVector<qint64> atms1;

        /** The force constant */
        SireUnits::Dimension::HarmonicBondConstant _k;

        /** The equilibrium distance for the restraint */
        SireUnits::Dimension::Length _r0;

        /** The equilibrium distance for the restraint at λ=0 */
        SireUnits::Dimension::Length _r1;
    };

    /** This class provides the information for a collection of Moving Harmonic Restraints
     *  restraints that can be added to a collection of molecues. Each
     *  restraint can act on a pair of particles or a pair of the
     *  centroids of two collections of particles.
     *  The restaints are spherically symmetric, and
     *  are simple harmonic potentials
     */
    class SIREMM_EXPORT MovingHarmonicRestraints
        : public SireBase::ConcreteProperty<MovingHarmonicRestraints, Restraints>
    {
        friend QDataStream & ::operator<<(QDataStream &, const SireMM::MovingHarmonicRestraints &);
        friend QDataStream & ::operator>>(QDataStream &, SireMM::MovingHarmonicRestraints &);

    public:
        MovingHarmonicRestraints();
        MovingHarmonicRestraints(const QString &name);

        MovingHarmonicRestraints(const MovingHarmonicRestraint &restraint);
        MovingHarmonicRestraints(const QList<MovingHarmonicRestraint> &restraints);

        MovingHarmonicRestraints(const QString &name,
                        const MovingHarmonicRestraint &restraint);

        MovingHarmonicRestraints(const QString &name,
                        const QList<MovingHarmonicRestraint> &restraints);

        MovingHarmonicRestraints(const MovingHarmonicRestraints &other);

        ~MovingHarmonicRestraints();

        MovingHarmonicRestraints &operator=(const MovingHarmonicRestraints &other);

        bool operator==(const MovingHarmonicRestraints &other) const;
        bool operator!=(const MovingHarmonicRestraints &other) const;

        static const char *typeName();
        const char *what() const;

        MovingHarmonicRestraints *clone() const;

        QString toString() const;

        bool isEmpty() const;
        bool isNull() const;

        int count() const;
        int size() const;
        int nRestraints() const;

        int nAtomRestraints() const;
        int nCentroidRestraints() const;

        bool hasAtomRestraints() const;
        bool hasCentroidRestraints() const;

        const MovingHarmonicRestraint &at(int i) const;
        const MovingHarmonicRestraint &operator[](int i) const;

        QList<MovingHarmonicRestraint> restraints() const;

        QList<MovingHarmonicRestraint> atomRestraints() const;
        QList<MovingHarmonicRestraint> centroidRestraints() const;

        void add(const MovingHarmonicRestraint &restraint);
        void add(const MovingHarmonicRestraints &restraints);

        MovingHarmonicRestraints &operator+=(const MovingHarmonicRestraint &restraint);
        MovingHarmonicRestraints &operator+=(const MovingHarmonicRestraints &restraints);

        MovingHarmonicRestraints operator+(const MovingHarmonicRestraint &restraint) const;
        MovingHarmonicRestraints operator+(const MovingHarmonicRestraints &restraints) const;

    private:
        /** The actual list of restraints*/
        QList<MovingHarmonicRestraint> r;
    };

}

Q_DECLARE_METATYPE(SireMM::MovingHarmonicRestraint)
Q_DECLARE_METATYPE(SireMM::MovingHarmonicRestraints)

SIRE_EXPOSE_CLASS(SireMM::MovingHarmonicRestraint)
SIRE_EXPOSE_CLASS(SireMM::MovingHarmonicRestraints)

SIRE_END_HEADER

#endif
