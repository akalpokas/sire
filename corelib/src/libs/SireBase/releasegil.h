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

#ifndef SIREBASE_RELEASEGIL_H
#define SIREBASE_RELEASEGIL_H

#include "sireglobal.h"
#include <memory>

SIRE_BEGIN_HEADER

namespace SireBase
{
    namespace detail
    {
        class ReleaseGILBase;
    }

    typedef std::shared_ptr<detail::ReleaseGILBase> GILHandle;

    SIREBASE_EXPORT GILHandle release_gil();

    namespace detail
    {
        /** This is the base class that will be sub-classed in the
         *  Python wrapper to provide the real code to release
         *  (and re-aquire) the GIL. We have to do this as
         *  the C++ libraries DO NOT link directly to Python. Only
         *  the wrappers link to Python.
         */
        class SIREBASE_EXPORT ReleaseGILBase
        {
            friend SIREBASE_EXPORT GILHandle SireBase::release_gil();

        public:
            ReleaseGILBase();
            virtual ~ReleaseGILBase();

        protected:
            static void registerReleaseGIL(ReleaseGILBase *handle);
            virtual GILHandle releaseGIL() const = 0;

        private:
            /** Pointer to the concrete class that is registered by
             *  the wrappers
             */
            static ReleaseGILBase *handle;
        };
    }
}

SIRE_END_HEADER

#endif
