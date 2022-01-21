// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Csch.pypp.hpp"

namespace bp = boost::python;

#include "SireStream/datastream.h"

#include "complexvalues.h"

#include "exp.h"

#include "expression.h"

#include "hyperbolicfuncs.h"

#include "identities.h"

#include "trigfuncs.h"

#include "hyperbolicfuncs.h"

SireCAS::Csch __copy__(const SireCAS::Csch &other){ return SireCAS::Csch(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_Csch_class(){

    { //::SireCAS::Csch
        typedef bp::class_< SireCAS::Csch, bp::bases< SireCAS::SingleFunc, SireCAS::ExBase > > Csch_exposer_t;
        Csch_exposer_t Csch_exposer = Csch_exposer_t( "Csch", "Hyperbolic cosecant", bp::init< >("Null constructor") );
        bp::scope Csch_scope( Csch_exposer );
        Csch_exposer.def( bp::init< SireCAS::Expression const & >(( bp::arg("ex") ), "Construct cos(expression)") );
        Csch_exposer.def( bp::init< SireCAS::Csch const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireCAS::Csch::evaluate
        
            typedef double ( ::SireCAS::Csch::*evaluate_function_type)( ::SireCAS::Values const & ) const;
            evaluate_function_type evaluate_function_value( &::SireCAS::Csch::evaluate );
            
            Csch_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , ( bp::arg("values") )
                , "Evaluate this function" );
        
        }
        { //::SireCAS::Csch::evaluate
        
            typedef ::SireMaths::Complex ( ::SireCAS::Csch::*evaluate_function_type)( ::SireCAS::ComplexValues const & ) const;
            evaluate_function_type evaluate_function_value( &::SireCAS::Csch::evaluate );
            
            Csch_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , ( bp::arg("values") )
                , "Complex evaluation" );
        
        }
        Csch_exposer.def( bp::self == bp::other< SireCAS::ExBase >() );
        { //::SireCAS::Csch::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireCAS::Csch::typeName );
            
            Csch_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireCAS::Csch::what
        
            typedef char const * ( ::SireCAS::Csch::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireCAS::Csch::what );
            
            Csch_exposer.def( 
                "what"
                , what_function_value
                , "" );
        
        }
        Csch_exposer.staticmethod( "typeName" );
        Csch_exposer.def( "__copy__", &__copy__);
        Csch_exposer.def( "__deepcopy__", &__copy__);
        Csch_exposer.def( "clone", &__copy__);
        Csch_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireCAS::Csch >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Csch_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireCAS::Csch >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Csch_exposer.def( "__setstate__", &__setstate__base64< ::SireCAS::Csch > );
        Csch_exposer.def( "__getstate__", &__getstate__base64< ::SireCAS::Csch > );
        Csch_exposer.def( "__str__", &__str__< ::SireCAS::Csch > );
        Csch_exposer.def( "__repr__", &__str__< ::SireCAS::Csch > );
        Csch_exposer.def( "__hash__", &::SireCAS::Csch::hash );
    }

}
