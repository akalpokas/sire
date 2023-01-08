// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "SymbolComplex.pypp.hpp"

namespace bp = boost::python;

#include "symbolcomplex.h"

#include "symbolcomplex.h"

SireCAS::SymbolComplex __copy__(const SireCAS::SymbolComplex &other){ return SireCAS::SymbolComplex(other); }

const char* pvt_get_name(const SireCAS::SymbolComplex&){ return "SireCAS::SymbolComplex";}

#include "Helpers/release_gil_policy.hpp"

void register_SymbolComplex_class(){

    { //::SireCAS::SymbolComplex
        typedef bp::class_< SireCAS::SymbolComplex > SymbolComplex_exposer_t;
        SymbolComplex_exposer_t SymbolComplex_exposer = SymbolComplex_exposer_t( "SymbolComplex", "Small class that holds a SymbolID number and an associated complex value", bp::init< SireCAS::SymbolID, SireMaths::Complex const & >(( bp::arg("id"), bp::arg("val") ), "") );
        bp::scope SymbolComplex_scope( SymbolComplex_exposer );
        SymbolComplex_exposer.def( bp::init< SireCAS::SymbolID, double >(( bp::arg("id"), bp::arg("val") ), "") );
        SymbolComplex_exposer.def( bp::init< SireCAS::SymbolValue const & >(( bp::arg("val") ), "") );
        { //::SireCAS::SymbolComplex::ID
        
            typedef ::SireCAS::SymbolID ( ::SireCAS::SymbolComplex::*ID_function_type)(  ) const;
            ID_function_type ID_function_value( &::SireCAS::SymbolComplex::ID );
            
            SymbolComplex_exposer.def( 
                "ID"
                , ID_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireCAS::SymbolComplex::value
        
            typedef ::SireMaths::Complex const & ( ::SireCAS::SymbolComplex::*value_function_type)(  ) const;
            value_function_type value_function_value( &::SireCAS::SymbolComplex::value );
            
            SymbolComplex_exposer.def( 
                "value"
                , value_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        SymbolComplex_exposer.def( "__copy__", &__copy__);
        SymbolComplex_exposer.def( "__deepcopy__", &__copy__);
        SymbolComplex_exposer.def( "clone", &__copy__);
        SymbolComplex_exposer.def( "__str__", &pvt_get_name);
        SymbolComplex_exposer.def( "__repr__", &pvt_get_name);
    }

}
