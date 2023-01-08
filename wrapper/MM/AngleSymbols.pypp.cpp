// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "AngleSymbols.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireFF/errors.h"

#include "SireMol/cgidx.h"

#include "SireMol/molecule.h"

#include "SireMol/mover.hpp"

#include "SireMol/partialmolecule.h"

#include "SireVol/coordgroup.h"

#include "internalparameters.h"

#include "sireglobal.h"

#include "tostring.h"

#include <algorithm>

#include "internalparameters.h"

SireMM::AngleSymbols __copy__(const SireMM::AngleSymbols &other){ return SireMM::AngleSymbols(other); }

const char* pvt_get_name(const SireMM::AngleSymbols&){ return "SireMM::AngleSymbols";}

#include "Helpers/release_gil_policy.hpp"

void register_AngleSymbols_class(){

    { //::SireMM::AngleSymbols
        typedef bp::class_< SireMM::AngleSymbols, bp::bases< SireMM::InternalSymbolsBase > > AngleSymbols_exposer_t;
        AngleSymbols_exposer_t AngleSymbols_exposer = AngleSymbols_exposer_t( "AngleSymbols", "This class holds the symbols required for the angle parameters.\nThese are functions of the angle 0-1-2 (theta) of three\natoms (where atom 1 is the central atom of the angle)", bp::init< >("") );
        bp::scope AngleSymbols_scope( AngleSymbols_exposer );
        { //::SireMM::AngleSymbols::theta
        
            typedef ::SireCAS::Symbol const & ( ::SireMM::AngleSymbols::*theta_function_type)(  ) const;
            theta_function_type theta_function_value( &::SireMM::AngleSymbols::theta );
            
            AngleSymbols_exposer.def( 
                "theta"
                , theta_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the symbols representing the angle (theta)" );
        
        }
        AngleSymbols_exposer.def( "__copy__", &__copy__);
        AngleSymbols_exposer.def( "__deepcopy__", &__copy__);
        AngleSymbols_exposer.def( "clone", &__copy__);
        AngleSymbols_exposer.def( "__str__", &pvt_get_name);
        AngleSymbols_exposer.def( "__repr__", &pvt_get_name);
    }

}
