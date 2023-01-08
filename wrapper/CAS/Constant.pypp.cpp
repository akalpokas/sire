// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "Constant.pypp.hpp"

namespace bp = boost::python;

#include "SireStream/datastream.h"

#include "complexvalues.h"

#include "constant.h"

#include "expression.h"

#include "expressions.h"

#include "functions.h"

#include "identities.h"

#include "integrationconstant.h"

#include "symbol.h"

#include "symbols.h"

#include "values.h"

#include "constant.h"

SireCAS::Constant __copy__(const SireCAS::Constant &other){ return SireCAS::Constant(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_Constant_class(){

    { //::SireCAS::Constant
        typedef bp::class_< SireCAS::Constant, bp::bases< SireCAS::ExBase > > Constant_exposer_t;
        Constant_exposer_t Constant_exposer = Constant_exposer_t( "Constant", "\nThis class represents a constant value (e.g. a number).\n\nAuthor: Christopher Woods\n", bp::init< >("Construct a constant") );
        bp::scope Constant_scope( Constant_exposer );
        Constant_exposer.def( bp::init< SireCAS::Constant const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireCAS::Constant::children
        
            typedef ::SireCAS::Expressions ( ::SireCAS::Constant::*children_function_type)(  ) const;
            children_function_type children_function_value( &::SireCAS::Constant::children );
            
            Constant_exposer.def( 
                "children"
                , children_function_value
                , bp::release_gil_policy()
                , "No children in a constant" );
        
        }
        { //::SireCAS::Constant::differentiate
        
            typedef ::SireCAS::Expression ( ::SireCAS::Constant::*differentiate_function_type)( ::SireCAS::Symbol const & ) const;
            differentiate_function_type differentiate_function_value( &::SireCAS::Constant::differentiate );
            
            Constant_exposer.def( 
                "differentiate"
                , differentiate_function_value
                , ( bp::arg("symbol") )
                , bp::release_gil_policy()
                , "Differential of a constant is zero" );
        
        }
        { //::SireCAS::Constant::evaluate
        
            typedef double ( ::SireCAS::Constant::*evaluate_function_type)( ::SireCAS::Values const & ) const;
            evaluate_function_type evaluate_function_value( &::SireCAS::Constant::evaluate );
            
            Constant_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , ( bp::arg("values") )
                , bp::release_gil_policy()
                , "Evaluation of a constant is 1" );
        
        }
        { //::SireCAS::Constant::evaluate
        
            typedef ::SireMaths::Complex ( ::SireCAS::Constant::*evaluate_function_type)( ::SireCAS::ComplexValues const & ) const;
            evaluate_function_type evaluate_function_value( &::SireCAS::Constant::evaluate );
            
            Constant_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , ( bp::arg("values") )
                , bp::release_gil_policy()
                , "Evaluation of a constant is 1" );
        
        }
        { //::SireCAS::Constant::expand
        
            typedef ::QList< SireCAS::Factor > ( ::SireCAS::Constant::*expand_function_type)( ::SireCAS::Symbol const & ) const;
            expand_function_type expand_function_value( &::SireCAS::Constant::expand );
            
            Constant_exposer.def( 
                "expand"
                , expand_function_value
                , ( bp::arg("symbol") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireCAS::Constant::functions
        
            typedef ::SireCAS::Functions ( ::SireCAS::Constant::*functions_function_type)(  ) const;
            functions_function_type functions_function_value( &::SireCAS::Constant::functions );
            
            Constant_exposer.def( 
                "functions"
                , functions_function_value
                , bp::release_gil_policy()
                , "No functions in a constant" );
        
        }
        { //::SireCAS::Constant::hash
        
            typedef ::uint ( ::SireCAS::Constant::*hash_function_type)(  ) const;
            hash_function_type hash_function_value( &::SireCAS::Constant::hash );
            
            Constant_exposer.def( 
                "hash"
                , hash_function_value
                , bp::release_gil_policy()
                , "Hash a constant" );
        
        }
        { //::SireCAS::Constant::integrate
        
            typedef ::SireCAS::Expression ( ::SireCAS::Constant::*integrate_function_type)( ::SireCAS::Symbol const & ) const;
            integrate_function_type integrate_function_value( &::SireCAS::Constant::integrate );
            
            Constant_exposer.def( 
                "integrate"
                , integrate_function_value
                , ( bp::arg("symbol") )
                , bp::release_gil_policy()
                , "Integral of a constant is = constantsymbol + C" );
        
        }
        Constant_exposer.def( bp::self == bp::other< SireCAS::ExBase >() );
        { //::SireCAS::Constant::substitute
        
            typedef ::SireCAS::Expression ( ::SireCAS::Constant::*substitute_function_type)( ::SireCAS::Identities const & ) const;
            substitute_function_type substitute_function_value( &::SireCAS::Constant::substitute );
            
            Constant_exposer.def( 
                "substitute"
                , substitute_function_value
                , ( bp::arg("identities") )
                , bp::release_gil_policy()
                , "Cant substitute into a constant" );
        
        }
        { //::SireCAS::Constant::symbols
        
            typedef ::SireCAS::Symbols ( ::SireCAS::Constant::*symbols_function_type)(  ) const;
            symbols_function_type symbols_function_value( &::SireCAS::Constant::symbols );
            
            Constant_exposer.def( 
                "symbols"
                , symbols_function_value
                , bp::release_gil_policy()
                , "No symbols in a constant" );
        
        }
        { //::SireCAS::Constant::toString
        
            typedef ::QString ( ::SireCAS::Constant::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireCAS::Constant::toString );
            
            Constant_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "Return a string representation of this constant (actually an empty string)" );
        
        }
        { //::SireCAS::Constant::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireCAS::Constant::typeName );
            
            Constant_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireCAS::Constant::what
        
            typedef char const * ( ::SireCAS::Constant::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireCAS::Constant::what );
            
            Constant_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        Constant_exposer.staticmethod( "typeName" );
        Constant_exposer.def( "__copy__", &__copy__);
        Constant_exposer.def( "__deepcopy__", &__copy__);
        Constant_exposer.def( "clone", &__copy__);
        Constant_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireCAS::Constant >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Constant_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireCAS::Constant >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Constant_exposer.def_pickle(sire_pickle_suite< ::SireCAS::Constant >());
        Constant_exposer.def( "__str__", &__str__< ::SireCAS::Constant > );
        Constant_exposer.def( "__repr__", &__str__< ::SireCAS::Constant > );
        Constant_exposer.def( "__hash__", &::SireCAS::Constant::hash );
    }

}
