// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "ExBase.pypp.hpp"

namespace bp = boost::python;

#include "SireCAS/errors.h"

#include "SireStream/datastream.h"

#include "complexvalues.h"

#include "exbase.h"

#include "expression.h"

#include "expressionbase.h"

#include "functions.h"

#include "identities.h"

#include "values.h"

#include "exbase.h"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_ExBase_class(){

    { //::SireCAS::ExBase
        typedef bp::class_< SireCAS::ExBase, boost::noncopyable > ExBase_exposer_t;
        ExBase_exposer_t ExBase_exposer = ExBase_exposer_t( "ExBase", "\nPure-virtual base class of all of the parts of mathematical expressions.\n\nThis class provides the atom of SireCAS. ExBase objects are combined together\nto form complete expressions. All constants, functions and symbols are derived\nfrom this object.\n\nThis class is an example of an implicitly shared, self-managed object, and\nis designed so that it can be held by SharedPolyPointer (indeed,\nExpressionBase is just a proxy for SharedPolyPointer<ExBase>).\n\nAuthor: Christopher Woods\n", bp::no_init );
        bp::scope ExBase_scope( ExBase_exposer );
        { //::SireCAS::ExBase::children
        
            typedef ::SireCAS::Expressions ( ::SireCAS::ExBase::*children_function_type)(  ) const;
            children_function_type children_function_value( &::SireCAS::ExBase::children );
            
            ExBase_exposer.def( 
                "children"
                , children_function_value
                , "Return the child expressions of this Expression" );
        
        }
        { //::SireCAS::ExBase::conjugate
        
            typedef ::SireCAS::Expression ( ::SireCAS::ExBase::*conjugate_function_type)(  ) const;
            conjugate_function_type conjugate_function_value( &::SireCAS::ExBase::conjugate );
            
            ExBase_exposer.def( 
                "conjugate"
                , conjugate_function_value
                , "Return the complex conjugate of this expression" );
        
        }
        { //::SireCAS::ExBase::differentiate
        
            typedef ::SireCAS::Expression ( ::SireCAS::ExBase::*differentiate_function_type)( ::SireCAS::Symbol const & ) const;
            differentiate_function_type differentiate_function_value( &::SireCAS::ExBase::differentiate );
            
            ExBase_exposer.def( 
                "differentiate"
                , differentiate_function_value
                , ( bp::arg("symbol") )
                , "Return an expression that the differential of this ExBase\nwith respect to symbol. Note an exception may\nbe thrown if this ExBase cannot be differentiated.\nThrow: SireCAS::unavailable_differential\n" );
        
        }
        { //::SireCAS::ExBase::evaluate
        
            typedef double ( ::SireCAS::ExBase::*evaluate_function_type)( ::SireCAS::Values const & ) const;
            evaluate_function_type evaluate_function_value( &::SireCAS::ExBase::evaluate );
            
            ExBase_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , ( bp::arg("values") )
                , "Evaluate this ExBase using values values. Any\nmissing symbols are assumed to equal zero.\n\nNote that an exception will be thrown if the result of the\nevaluation of this, or one of its children, is complex.\n\nThrow: SireMaths::domain_error\n" );
        
        }
        { //::SireCAS::ExBase::evaluate
        
            typedef ::SireMaths::Complex ( ::SireCAS::ExBase::*evaluate_function_type)( ::SireCAS::ComplexValues const & ) const;
            evaluate_function_type evaluate_function_value( &::SireCAS::ExBase::evaluate );
            
            ExBase_exposer.def( 
                "evaluate"
                , evaluate_function_value
                , ( bp::arg("values") )
                , "Evaluate this ExBase using the complex values values.\nAny missing symbols are assumed to equal zero." );
        
        }
        { //::SireCAS::ExBase::expand
        
            typedef ::QList< SireCAS::Factor > ( ::SireCAS::ExBase::*expand_function_type)( ::SireCAS::Symbol const & ) const;
            expand_function_type expand_function_value( &::SireCAS::ExBase::expand );
            
            ExBase_exposer.def( 
                "expand"
                , expand_function_value
                , ( bp::arg("symbol") )
                , "Rearrange this expression into the form\nm x^i + n x^j + ... + constant\nand return the factors and powers of x\n\nThrow: SireCAS::rearrangement_error\n" );
        
        }
        { //::SireCAS::ExBase::functions
        
            typedef ::SireCAS::Functions ( ::SireCAS::ExBase::*functions_function_type)(  ) const;
            functions_function_type functions_function_value( &::SireCAS::ExBase::functions );
            
            ExBase_exposer.def( 
                "functions"
                , functions_function_value
                , "Return the set of Functions that appear in this ExBase" );
        
        }
        { //::SireCAS::ExBase::hash
        
            typedef ::uint ( ::SireCAS::ExBase::*hash_function_type)(  ) const;
            hash_function_type hash_function_value( &::SireCAS::ExBase::hash );
            
            ExBase_exposer.def( 
                "hash"
                , hash_function_value
                , "Return a hash of this object - return a combination of the\nidentifying magic for the class and a hash for its contents." );
        
        }
        { //::SireCAS::ExBase::integrate
        
            typedef ::SireCAS::Expression ( ::SireCAS::ExBase::*integrate_function_type)( ::SireCAS::Symbol const & ) const;
            integrate_function_type integrate_function_value( &::SireCAS::ExBase::integrate );
            
            ExBase_exposer.def( 
                "integrate"
                , integrate_function_value
                , ( bp::arg("symbol") )
                , "Return the indefinite integral of this ExBase with respect to\nsymbol. This is not guaranteed to work() and will return an\nexpression of the form Sum( integral(exbase) + integral_constant ).\nIf it doesnt work then an exception will be throw.\nThrow: SireCAS::unavailable_integral\n" );
        
        }
        { //::SireCAS::ExBase::isComplex
        
            typedef bool ( ::SireCAS::ExBase::*isComplex_function_type)(  ) const;
            isComplex_function_type isComplex_function_value( &::SireCAS::ExBase::isComplex );
            
            ExBase_exposer.def( 
                "isComplex"
                , isComplex_function_value
                , "Return whether or not this expression contains any complex (imaginary)\nparts" );
        
        }
        { //::SireCAS::ExBase::isCompound
        
            typedef bool ( ::SireCAS::ExBase::*isCompound_function_type)(  ) const;
            isCompound_function_type isCompound_function_value( &::SireCAS::ExBase::isCompound );
            
            ExBase_exposer.def( 
                "isCompound"
                , isCompound_function_value
                , "Return whether or not this is a compound expression, and thus as such\nrequires brackets placed around it when it is printed. Examples include\nSum, Product and Power. For most other functions it is safe to leave\nthis as false." );
        
        }
        { //::SireCAS::ExBase::isConstant
        
            typedef bool ( ::SireCAS::ExBase::*isConstant_function_type)(  ) const;
            isConstant_function_type isConstant_function_value( &::SireCAS::ExBase::isConstant );
            
            ExBase_exposer.def( 
                "isConstant"
                , isConstant_function_value
                , "Return whether or not this is a constant expression (does not\ndepend on any symbols)" );
        
        }
        { //::SireCAS::ExBase::isFunction
        
            typedef bool ( ::SireCAS::ExBase::*isFunction_function_type)( ::SireCAS::Symbol const & ) const;
            isFunction_function_type isFunction_function_value( &::SireCAS::ExBase::isFunction );
            
            ExBase_exposer.def( 
                "isFunction"
                , isFunction_function_value
                , ( bp::arg("arg0") )
                , "Return whether or not this is a function of the passed Symbol" );
        
        }
        ExBase_exposer.def( bp::self != bp::self );
        ExBase_exposer.def( -bp::self );
        { //::SireCAS::ExBase::operator=
        
            typedef ::SireCAS::ExBase & ( ::SireCAS::ExBase::*assign_function_type)( ::SireCAS::ExBase const & ) ;
            assign_function_type assign_function_value( &::SireCAS::ExBase::operator= );
            
            ExBase_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("arg0") )
                , bp::return_self< >()
                , "Assignment operator" );
        
        }
        ExBase_exposer.def( bp::self == bp::self );
        { //::SireCAS::ExBase::series
        
            typedef ::SireCAS::Expression ( ::SireCAS::ExBase::*series_function_type)( ::SireCAS::Symbol const &,int ) const;
            series_function_type series_function_value( &::SireCAS::ExBase::series );
            
            ExBase_exposer.def( 
                "series"
                , series_function_value
                , ( bp::arg("symbol"), bp::arg("n") )
                , "Return a series expansion of this expression with respect to\nsymbol, if possible, to order\nn. This is not guaranteed to work, and will return this expression\nunchanged if it doesnt work. If it is expanded, then a series\nwill be returned, together with an estimate of the error (e.g. O(x^2))" );
        
        }
        { //::SireCAS::ExBase::simplify
        
            typedef ::SireCAS::Expression ( ::SireCAS::ExBase::*simplify_function_type)( int ) const;
            simplify_function_type simplify_function_value( &::SireCAS::ExBase::simplify );
            
            ExBase_exposer.def( 
                "simplify"
                , simplify_function_value
                , ( bp::arg("options")=(int)(0) )
                , "Try to simplify this expression. This will try to use known mathematical\nidentities to convert complex expressions down to more simple ones.\nIf SireCAS::UNSAFE_COMPLEX_SIMPLIFICATIONS is true, then identities\nthat are not safe for complex math are used, e.g. z = sin(arcsin(z))." );
        
        }
        { //::SireCAS::ExBase::substitute
        
            typedef ::SireCAS::Expression ( ::SireCAS::ExBase::*substitute_function_type)( ::SireCAS::Identities const & ) const;
            substitute_function_type substitute_function_value( &::SireCAS::ExBase::substitute );
            
            ExBase_exposer.def( 
                "substitute"
                , substitute_function_value
                , ( bp::arg("identities") )
                , "Return an expression that has the identities in identities\nsubstituted into this expression" );
        
        }
        { //::SireCAS::ExBase::symbols
        
            typedef ::SireCAS::Symbols ( ::SireCAS::ExBase::*symbols_function_type)(  ) const;
            symbols_function_type symbols_function_value( &::SireCAS::ExBase::symbols );
            
            ExBase_exposer.def( 
                "symbols"
                , symbols_function_value
                , "Return the set of Symbols that appear in this ExBase" );
        
        }
        { //::SireCAS::ExBase::toOpenMMString
        
            typedef ::QString ( ::SireCAS::ExBase::*toOpenMMString_function_type)(  ) const;
            toOpenMMString_function_type toOpenMMString_function_value( &::SireCAS::ExBase::toOpenMMString );
            
            ExBase_exposer.def( 
                "toOpenMMString"
                , toOpenMMString_function_value
                , "Return a string representation of this object in the OpenMM syntax" );
        
        }
        { //::SireCAS::ExBase::toString
        
            typedef ::QString ( ::SireCAS::ExBase::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireCAS::ExBase::toString );
            
            ExBase_exposer.def( 
                "toString"
                , toString_function_value
                , "Return a string representation of this object" );
        
        }
        { //::SireCAS::ExBase::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireCAS::ExBase::typeName );
            
            ExBase_exposer.def( 
                "typeName"
                , typeName_function_value
                , "Return the name of this class type" );
        
        }
        { //::SireCAS::ExBase::what
        
            typedef char const * ( ::SireCAS::ExBase::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireCAS::ExBase::what );
            
            ExBase_exposer.def( 
                "what"
                , what_function_value
                , "Return the name of the type of this ExBase object" );
        
        }
        ExBase_exposer.staticmethod( "typeName" );
        ExBase_exposer.def( bp::self * bp::self );
        ExBase_exposer.def( bp::self * bp::other< SireCAS::Expression >() );
        ExBase_exposer.def( bp::self * bp::other< double >() );
        ExBase_exposer.def( bp::other< double >() * bp::self );
        ExBase_exposer.def( bp::self * bp::other< SireMaths::Complex >() );
        ExBase_exposer.def( bp::self + bp::self );
        ExBase_exposer.def( bp::self + bp::other< SireCAS::Expression >() );
        ExBase_exposer.def( bp::self + bp::other< double >() );
        ExBase_exposer.def( bp::other< double >() + bp::self );
        ExBase_exposer.def( bp::self + bp::other< SireMaths::Complex >() );
        ExBase_exposer.def( bp::self - bp::self );
        ExBase_exposer.def( bp::self - bp::other< SireCAS::Expression >() );
        ExBase_exposer.def( bp::self - bp::other< double >() );
        ExBase_exposer.def( bp::other< double >() - bp::self );
        ExBase_exposer.def( bp::self - bp::other< SireMaths::Complex >() );
        ExBase_exposer.def( bp::self / bp::self );
        ExBase_exposer.def( bp::self / bp::other< SireCAS::Expression >() );
        ExBase_exposer.def( bp::self / bp::other< double >() );
        ExBase_exposer.def( bp::other< double >() / bp::self );
        ExBase_exposer.def( bp::self / bp::other< SireMaths::Complex >() );
        ExBase_exposer.def( self + self );
        ExBase_exposer.def( self - self );
        ExBase_exposer.def( self * self );
        ExBase_exposer.def( self / self );
        ExBase_exposer.def( other<double>() + self );
        ExBase_exposer.def( other<double>() - self );
        ExBase_exposer.def( other<double>() * self );
        ExBase_exposer.def( other<double>() / self );
        ExBase_exposer.def( self + other<double>() );
        ExBase_exposer.def( self - other<double>() );
        ExBase_exposer.def( self * other<double>() );
        ExBase_exposer.def( self / other<double>() );
        ExBase_exposer.def( other<SireMaths::Complex>() + self );
        ExBase_exposer.def( other<SireMaths::Complex>() - self );
        ExBase_exposer.def( other<SireMaths::Complex>() * self );
        ExBase_exposer.def( other<SireMaths::Complex>() / self );
        ExBase_exposer.def( self + other<SireMaths::Complex>() );
        ExBase_exposer.def( self - other<SireMaths::Complex>() );
        ExBase_exposer.def( self * other<SireMaths::Complex>() );
        ExBase_exposer.def( self / other<SireMaths::Complex>() );
        ExBase_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireCAS::ExBase >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ExBase_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireCAS::ExBase >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ExBase_exposer.def( "__setstate__", &__setstate__base64< ::SireCAS::ExBase > );
        ExBase_exposer.def( "__getstate__", &__getstate__base64< ::SireCAS::ExBase > );
        ExBase_exposer.def( "__str__", &__str__< ::SireCAS::ExBase > );
        ExBase_exposer.def( "__repr__", &__str__< ::SireCAS::ExBase > );
        ExBase_exposer.def( "__hash__", &::SireCAS::ExBase::hash );
    }

}
