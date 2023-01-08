// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "ImproperSymbols.pypp.hpp"

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

SireMM::ImproperSymbols __copy__(const SireMM::ImproperSymbols &other){ return SireMM::ImproperSymbols(other); }

const char* pvt_get_name(const SireMM::ImproperSymbols&){ return "SireMM::ImproperSymbols";}

#include "Helpers/release_gil_policy.hpp"

void register_ImproperSymbols_class(){

    { //::SireMM::ImproperSymbols
        typedef bp::class_< SireMM::ImproperSymbols, bp::bases< SireMM::InternalSymbolsBase > > ImproperSymbols_exposer_t;
        ImproperSymbols_exposer_t ImproperSymbols_exposer = ImproperSymbols_exposer_t( "ImproperSymbols", "This class holds the symbols required for the improper parameters.\nThese are functions of the dihedral (0-1-2-3) (phi) and also\nof the improper angle made between the bond 0-1 and the plane\nformed by atoms 1-2-3 (theta)", bp::init< >("") );
        bp::scope ImproperSymbols_scope( ImproperSymbols_exposer );
        { //::SireMM::ImproperSymbols::phi
        
            typedef ::SireCAS::Symbol const & ( ::SireMM::ImproperSymbols::*phi_function_type)(  ) const;
            phi_function_type phi_function_value( &::SireMM::ImproperSymbols::phi );
            
            ImproperSymbols_exposer.def( 
                "phi"
                , phi_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the symbol representing the torsion 0-1-2-3" );
        
        }
        { //::SireMM::ImproperSymbols::theta
        
            typedef ::SireCAS::Symbol const & ( ::SireMM::ImproperSymbols::*theta_function_type)(  ) const;
            theta_function_type theta_function_value( &::SireMM::ImproperSymbols::theta );
            
            ImproperSymbols_exposer.def( 
                "theta"
                , theta_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the symbol representing the angle between the improper\nand the plane formed by atoms 1-3" );
        
        }
        ImproperSymbols_exposer.def( "__copy__", &__copy__);
        ImproperSymbols_exposer.def( "__deepcopy__", &__copy__);
        ImproperSymbols_exposer.def( "clone", &__copy__);
        ImproperSymbols_exposer.def( "__str__", &pvt_get_name);
        ImproperSymbols_exposer.def( "__repr__", &pvt_get_name);
    }

}
