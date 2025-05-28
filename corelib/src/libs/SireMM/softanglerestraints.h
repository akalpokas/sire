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

#ifndef SIREMM_SOFTANGLERESTRAINTS_H
#define SIREMM_SOFTANGLERESTRAINTS_H


#include "restraints.h"
#include "SireUnits/dimensions.h"
#include "SireUnits/generalunit.h"

SIRE_BEGIN_HEADER

namespace SireMM
{
    class SoftAngleRestraint;
    class SoftAngleRestraints;
}

SIREMM_EXPORT QDataStream &operator<<(QDataStream &, const SireMM::SoftAngleRestraint &);
SIREMM_EXPORT QDataStream &operator>>(QDataStream &, SireMM::SoftAngleRestraint &);

SIREMM_EXPORT QDataStream &operator<<(QDataStream &, const SireMM::SoftAngleRestraints &);
SIREMM_EXPORT QDataStream &operator>>(QDataStream &, SireMM::SoftAngleRestraints &);

namespace SireMM
{

    /** This class represents a single angle restraint between any three
     *  atoms in a system softened by a plateau potential.
     * @author Audrius Kalpokas
     */
    class SIREMM_EXPORT SoftAngleRestraint
        : public SireBase::ConcreteProperty<SoftAngleRestraint, SireBase::Property>
    {

        friend SIREMM_EXPORT QDataStream & ::operator<<(QDataStream &, const SireMM::SoftAngleRestraint &);
        friend SIREMM_EXPORT QDataStream & ::operator>>(QDataStream &, SireMM::SoftAngleRestraint &);

    public:
        SoftAngleRestraint();
        SoftAngleRestraint(const QList<qint64> &atoms,
                       const SireUnits::Dimension::Angle &theta0,
                       const SireUnits::Dimension::HarmonicAngleConstant &ktheta,
                       const SireUnits::Dimension::MolarEnergy &pe);


        SoftAngleRestraint(const SoftAngleRestraint &other);

        ~SoftAngleRestraint();

        SoftAngleRestraint &operator=(const SoftAngleRestraint &other);

        bool operator==(const SoftAngleRestraint &other) const;
        bool operator!=(const SoftAngleRestraint &other) const;

        SoftAngleRestraints operator+(const SoftAngleRestraint &other) const;
        SoftAngleRestraints operator+(const SoftAngleRestraints &other) const;

        static const char *typeName();
        const char *what() const;

        SoftAngleRestraint *clone() const;

        QString toString() const;

        bool isNull() const;

        QVector<qint64> atoms() const;

        SireUnits::Dimension::Angle theta0() const;
        SireUnits::Dimension::HarmonicAngleConstant ktheta() const;
        SireUnits::Dimension::MolarEnergy pe() const;

    private:
        /** Atoms involved in the angle restraint */
        QVector<qint64> atms;

        /** Equilibrium angle */
        SireUnits::Dimension::Angle _theta0;

        /** Harmonic angle constant */
        SireUnits::Dimension::HarmonicAngleConstant _ktheta;

        /** Plateau energy of the restraint */
        SireUnits::Dimension::MolarEnergy _pe;
    };

    /** This class represents a collection of angle restraints */
    class SIREMM_EXPORT SoftAngleRestraints
        : public SireBase::ConcreteProperty<SoftAngleRestraints, Restraints>
    {
        friend SIREMM_EXPORT QDataStream & ::operator<<(QDataStream &, const SireMM::SoftAngleRestraints &);
        friend SIREMM_EXPORT QDataStream & ::operator>>(QDataStream &, SireMM::SoftAngleRestraints &);

    public:
        SoftAngleRestraints();
        SoftAngleRestraints(const QString &name);

        SoftAngleRestraints(const SoftAngleRestraint &restraint);
        SoftAngleRestraints(const QList<SoftAngleRestraint> &restraints);

        SoftAngleRestraints(const QString &name,
                        const SoftAngleRestraint &restraint);

        SoftAngleRestraints(const QString &name,
                        const QList<SoftAngleRestraint> &restraints);

        SoftAngleRestraints(const SoftAngleRestraints &other);

        ~SoftAngleRestraints();

        SoftAngleRestraints &operator=(const SoftAngleRestraints &other);

        bool operator==(const SoftAngleRestraints &other) const;
        bool operator!=(const SoftAngleRestraints &other) const;

        static const char *typeName();
        const char *what() const;

        SoftAngleRestraints *clone() const;

        QString toString() const;

        bool isEmpty() const;
        bool isNull() const;

        int count() const;
        int size() const;
        int nRestraints() const;

        const SoftAngleRestraint &at(int i) const;
        const SoftAngleRestraint &operator[](int i) const;

        QList<SoftAngleRestraint> restraints() const;

        void add(const SoftAngleRestraint &restraint);
        void add(const SoftAngleRestraints &restraints);

        SoftAngleRestraints &operator+=(const SoftAngleRestraint &restraint);
        SoftAngleRestraints &operator+=(const SoftAngleRestraints &restraints);

        SoftAngleRestraints operator+(const SoftAngleRestraint &restraint) const;
        SoftAngleRestraints operator+(const SoftAngleRestraints &restraints) const;

    private:
        /** List of restraints */
        QList<SoftAngleRestraint> r;
    };
}

Q_DECLARE_METATYPE(SireMM::SoftAngleRestraint)
Q_DECLARE_METATYPE(SireMM::SoftAngleRestraints)

SIRE_EXPOSE_CLASS(SireMM::SoftAngleRestraint)
SIRE_EXPOSE_CLASS(SireMM::SoftAngleRestraints)
SIRE_END_HEADER

#endif
