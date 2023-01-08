// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "StretchBendTorsionSymbols.pypp.hpp"

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

SireMM::StretchBendTorsionSymbols __copy__(const SireMM::StretchBendTorsionSymbols &other){ return SireMM::StretchBendTorsionSymbols(other); }

const char* pvt_get_name(const SireMM::StretchBendTorsionSymbols&){ return "SireMM::StretchBendTorsionSymbols";}

#include "Helpers/release_gil_policy.hpp"

void register_StretchBendTorsionSymbols_class(){

    { //::SireMM::StretchBendTorsionSymbols
        typedef bp::class_< SireMM::StretchBendTorsionSymbols, bp::bases< SireMM::InternalSymbolsBase > > StretchBendTorsionSymbols_exposer_t;
        StretchBendTorsionSymbols_exposer_t StretchBendTorsionSymbols_exposer = StretchBendTorsionSymbols_exposer_t( "StretchBendTorsionSymbols", "This class holds the symbols required for the stretch-bend-torsion parameters\nThese are functions over four atoms, 0, 1, 2, 3 and 4, of the dihedral\nformed over the four atoms 0-1-2-3 (phi), the angles 0-1-2 (theta012)\nand 3-2-1 (theta321) and the distances between atoms, 0-1 (r01),\n1-2 (r12), 3-2 (r32) and 0-3 (r03)\n", bp::init< >("") );
        bp::scope StretchBendTorsionSymbols_scope( StretchBendTorsionSymbols_exposer );
        { //::SireMM::StretchBendTorsionSymbols::phi
        
            typedef ::SireCAS::Symbol const & ( ::SireMM::StretchBendTorsionSymbols::*phi_function_type)(  ) const;
            phi_function_type phi_function_value( &::SireMM::StretchBendTorsionSymbols::phi );
            
            StretchBendTorsionSymbols_exposer.def( 
                "phi"
                , phi_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the symbol representing the torsion, phi" );
        
        }
        { //::SireMM::StretchBendTorsionSymbols::r01
        
            typedef ::SireCAS::Symbol const & ( ::SireMM::StretchBendTorsionSymbols::*r01_function_type)(  ) const;
            r01_function_type r01_function_value( &::SireMM::StretchBendTorsionSymbols::r01 );
            
            StretchBendTorsionSymbols_exposer.def( 
                "r01"
                , r01_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the symbol representing the bond between atoms 0-1, r_\n{01}" );
        
        }
        { //::SireMM::StretchBendTorsionSymbols::r03
        
            typedef ::SireCAS::Symbol const & ( ::SireMM::StretchBendTorsionSymbols::*r03_function_type)(  ) const;
            r03_function_type r03_function_value( &::SireMM::StretchBendTorsionSymbols::r03 );
            
            StretchBendTorsionSymbols_exposer.def( 
                "r03"
                , r03_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the symbol representing the distance from atom 0 to 3, r_\n{03}" );
        
        }
        { //::SireMM::StretchBendTorsionSymbols::r12
        
            typedef ::SireCAS::Symbol const & ( ::SireMM::StretchBendTorsionSymbols::*r12_function_type)(  ) const;
            r12_function_type r12_function_value( &::SireMM::StretchBendTorsionSymbols::r12 );
            
            StretchBendTorsionSymbols_exposer.def( 
                "r12"
                , r12_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the symbol representing the bond between atoms 1-2, r_\n{12}" );
        
        }
        { //::SireMM::StretchBendTorsionSymbols::r32
        
            typedef ::SireCAS::Symbol const & ( ::SireMM::StretchBendTorsionSymbols::*r32_function_type)(  ) const;
            r32_function_type r32_function_value( &::SireMM::StretchBendTorsionSymbols::r32 );
            
            StretchBendTorsionSymbols_exposer.def( 
                "r32"
                , r32_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the symbol representing the bond between atoms 3-2, r_\n{32}" );
        
        }
        { //::SireMM::StretchBendTorsionSymbols::theta012
        
            typedef ::SireCAS::Symbol const & ( ::SireMM::StretchBendTorsionSymbols::*theta012_function_type)(  ) const;
            theta012_function_type theta012_function_value( &::SireMM::StretchBendTorsionSymbols::theta012 );
            
            StretchBendTorsionSymbols_exposer.def( 
                "theta012"
                , theta012_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the symbol representing the angle between atoms 0-1-2, theta_\n{012}" );
        
        }
        { //::SireMM::StretchBendTorsionSymbols::theta321
        
            typedef ::SireCAS::Symbol const & ( ::SireMM::StretchBendTorsionSymbols::*theta321_function_type)(  ) const;
            theta321_function_type theta321_function_value( &::SireMM::StretchBendTorsionSymbols::theta321 );
            
            StretchBendTorsionSymbols_exposer.def( 
                "theta321"
                , theta321_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the symbol representing the angle between atoms 3-2-1, theta_\n{321}" );
        
        }
        StretchBendTorsionSymbols_exposer.def( "__copy__", &__copy__);
        StretchBendTorsionSymbols_exposer.def( "__deepcopy__", &__copy__);
        StretchBendTorsionSymbols_exposer.def( "clone", &__copy__);
        StretchBendTorsionSymbols_exposer.def( "__str__", &pvt_get_name);
        StretchBendTorsionSymbols_exposer.def( "__repr__", &pvt_get_name);
    }

}
