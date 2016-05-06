// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "ArcTanh.pypp.hpp"

namespace bp = boost::python;

#include "SireMaths/errors.h"

#include "SireStream/datastream.h"

#include "complexvalues.h"

#include "exp.h"

#include "expression.h"

#include "hyperbolicfuncs.h"

#include "identities.h"

#include "invhyperbolicfuncs.h"

#include "invtrigfuncs.h"

#include "trigfuncs.h"

#include "invhyperbolicfuncs.h"

SireCAS::ArcTanh __copy__(const SireCAS::ArcTanh &other){ return SireCAS::ArcTanh(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_ArcTanh_class(){

    { //::SireCAS::ArcTanh
        typedef bp::class_< SireCAS::ArcTanh, bp::bases< SireCAS::SingleFunc, SireCAS::ExBase > > ArcTanh_exposer_t;
        ArcTanh_exposer_t ArcTanh_exposer = ArcTanh_exposer_t( "ArcTanh", bp::init< >() );
        bp::scope ArcTanh_scope( ArcTanh_exposer );
        ArcTanh_exposer.def( bp::init< SireCAS::Expression const & >(( bp::arg("ex") )) );
        ArcTanh_exposer.def( bp::init< SireCAS::ArcTanh const & >(( bp::arg("other") )) );
        { //::SireCAS::ArcTanh::evaluate
        
            typedef double ( ::SireCAS::ArcTanh::*evaluate_function_type)( ::SireCAS::Values const & ) const;
            evaluate_function_type evaluate_function_value( &::SireCAS::ArcTanh::evaluate );
            
            ArcTanh_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , ( bp::arg("values") ) );
        
        }
        { //::SireCAS::ArcTanh::evaluate
        
            typedef ::SireMaths::Complex ( ::SireCAS::ArcTanh::*evaluate_function_type)( ::SireCAS::ComplexValues const & ) const;
            evaluate_function_type evaluate_function_value( &::SireCAS::ArcTanh::evaluate );
            
            ArcTanh_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , ( bp::arg("values") ) );
        
        }
        ArcTanh_exposer.def( bp::self == bp::other< SireCAS::ExBase >() );
        { //::SireCAS::ArcTanh::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireCAS::ArcTanh::typeName );
            
            ArcTanh_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireCAS::ArcTanh::what
        
            typedef char const * ( ::SireCAS::ArcTanh::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireCAS::ArcTanh::what );
            
            ArcTanh_exposer.def( 
                "what"
                , what_function_value );
        
        }
        ArcTanh_exposer.staticmethod( "typeName" );
        ArcTanh_exposer.def( "__copy__", &__copy__);
        ArcTanh_exposer.def( "__deepcopy__", &__copy__);
        ArcTanh_exposer.def( "clone", &__copy__);
        ArcTanh_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireCAS::ArcTanh >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ArcTanh_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireCAS::ArcTanh >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ArcTanh_exposer.def( "__str__", &__str__< ::SireCAS::ArcTanh > );
        ArcTanh_exposer.def( "__repr__", &__str__< ::SireCAS::ArcTanh > );
        ArcTanh_exposer.def( "__hash__", &::SireCAS::ArcTanh::hash );
    }

}
