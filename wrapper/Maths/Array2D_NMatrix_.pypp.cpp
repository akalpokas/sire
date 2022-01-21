// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Array2D_NMatrix_.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/array2d.hpp"

#include "SireError/errors.h"

#include "SireMaths/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "n4matrix.h"

#include "nmatrix.h"

#include "nvector.h"

#include "n4matrix.h"

SireBase::Array2D<SireMaths::NMatrix> __copy__(const SireBase::Array2D<SireMaths::NMatrix> &other){ return SireBase::Array2D<SireMaths::NMatrix>(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_Array2D_NMatrix__class(){

    { //::SireBase::Array2D< SireMaths::NMatrix >
        typedef bp::class_< SireBase::Array2D< SireMaths::NMatrix >, bp::bases< SireBase::Array2DBase > > Array2D_NMatrix__exposer_t;
        Array2D_NMatrix__exposer_t Array2D_NMatrix__exposer = Array2D_NMatrix__exposer_t( "Array2D_NMatrix_", "", bp::init< >("") );
        bp::scope Array2D_NMatrix__scope( Array2D_NMatrix__exposer );
        Array2D_NMatrix__exposer.def( bp::init< int, int >(( bp::arg("nrows"), bp::arg("ncolumns") ), "") );
        Array2D_NMatrix__exposer.def( bp::init< int, int, SireMaths::NMatrix const & >(( bp::arg("nrows"), bp::arg("ncolumns"), bp::arg("default_value") ), "") );
        Array2D_NMatrix__exposer.def( bp::init< SireBase::Array2D< SireMaths::NMatrix > const & >(( bp::arg("other") ), "") );
        { //::SireBase::Array2D< SireMaths::NMatrix >::at
        
            typedef SireBase::Array2D< SireMaths::NMatrix > exported_class_t;
            typedef ::SireMaths::NMatrix const & ( ::SireBase::Array2D< SireMaths::NMatrix >::*at_function_type)( int,int ) const;
            at_function_type at_function_value( &::SireBase::Array2D< SireMaths::NMatrix >::at );
            
            Array2D_NMatrix__exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("i"), bp::arg("j") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireBase::Array2D< SireMaths::NMatrix >::get
        
            typedef SireBase::Array2D< SireMaths::NMatrix > exported_class_t;
            typedef ::SireMaths::NMatrix const & ( ::SireBase::Array2D< SireMaths::NMatrix >::*get_function_type)( int,int ) const;
            get_function_type get_function_value( &::SireBase::Array2D< SireMaths::NMatrix >::get );
            
            Array2D_NMatrix__exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("i"), bp::arg("j") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        Array2D_NMatrix__exposer.def( bp::self != bp::self );
        { //::SireBase::Array2D< SireMaths::NMatrix >::operator()
        
            typedef SireBase::Array2D< SireMaths::NMatrix > exported_class_t;
            typedef ::SireMaths::NMatrix const & ( ::SireBase::Array2D< SireMaths::NMatrix >::*__call___function_type)( int,int ) const;
            __call___function_type __call___function_value( &::SireBase::Array2D< SireMaths::NMatrix >::operator() );
            
            Array2D_NMatrix__exposer.def( 
                "__call__"
                , __call___function_value
                , ( bp::arg("i"), bp::arg("j") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireBase::Array2D< SireMaths::NMatrix >::operator=
        
            typedef SireBase::Array2D< SireMaths::NMatrix > exported_class_t;
            typedef ::SireBase::Array2D< SireMaths::NMatrix > & ( ::SireBase::Array2D< SireMaths::NMatrix >::*assign_function_type)( ::SireBase::Array2D< SireMaths::NMatrix > const & ) ;
            assign_function_type assign_function_value( &::SireBase::Array2D< SireMaths::NMatrix >::operator= );
            
            Array2D_NMatrix__exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        Array2D_NMatrix__exposer.def( bp::self == bp::self );
        { //::SireBase::Array2D< SireMaths::NMatrix >::redimension
        
            typedef SireBase::Array2D< SireMaths::NMatrix > exported_class_t;
            typedef void ( ::SireBase::Array2D< SireMaths::NMatrix >::*redimension_function_type)( int,int ) ;
            redimension_function_type redimension_function_value( &::SireBase::Array2D< SireMaths::NMatrix >::redimension );
            
            Array2D_NMatrix__exposer.def( 
                "redimension"
                , redimension_function_value
                , ( bp::arg("nrows"), bp::arg("ncolumns") )
                , "" );
        
        }
        { //::SireBase::Array2D< SireMaths::NMatrix >::set
        
            typedef SireBase::Array2D< SireMaths::NMatrix > exported_class_t;
            typedef void ( ::SireBase::Array2D< SireMaths::NMatrix >::*set_function_type)( int,int,::SireMaths::NMatrix const & ) ;
            set_function_type set_function_value( &::SireBase::Array2D< SireMaths::NMatrix >::set );
            
            Array2D_NMatrix__exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("i"), bp::arg("j"), bp::arg("value") )
                , "" );
        
        }
        { //::SireBase::Array2D< SireMaths::NMatrix >::setAll
        
            typedef SireBase::Array2D< SireMaths::NMatrix > exported_class_t;
            typedef void ( ::SireBase::Array2D< SireMaths::NMatrix >::*setAll_function_type)( ::SireMaths::NMatrix const & ) ;
            setAll_function_type setAll_function_value( &::SireBase::Array2D< SireMaths::NMatrix >::setAll );
            
            Array2D_NMatrix__exposer.def( 
                "setAll"
                , setAll_function_value
                , ( bp::arg("value") )
                , "" );
        
        }
        { //::SireBase::Array2D< SireMaths::NMatrix >::toString
        
            typedef SireBase::Array2D< SireMaths::NMatrix > exported_class_t;
            typedef ::QString ( ::SireBase::Array2D< SireMaths::NMatrix >::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireBase::Array2D< SireMaths::NMatrix >::toString );
            
            Array2D_NMatrix__exposer.def( 
                "toString"
                , toString_function_value
                , "" );
        
        }
        { //::SireBase::Array2D< SireMaths::NMatrix >::transpose
        
            typedef SireBase::Array2D< SireMaths::NMatrix > exported_class_t;
            typedef ::SireBase::Array2D< SireMaths::NMatrix > ( ::SireBase::Array2D< SireMaths::NMatrix >::*transpose_function_type)(  ) const;
            transpose_function_type transpose_function_value( &::SireBase::Array2D< SireMaths::NMatrix >::transpose );
            
            Array2D_NMatrix__exposer.def( 
                "transpose"
                , transpose_function_value
                , "" );
        
        }
        Array2D_NMatrix__exposer.def( "__copy__", &__copy__);
        Array2D_NMatrix__exposer.def( "__deepcopy__", &__copy__);
        Array2D_NMatrix__exposer.def( "clone", &__copy__);
        Array2D_NMatrix__exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireBase::Array2D<SireMaths::NMatrix> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Array2D_NMatrix__exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireBase::Array2D<SireMaths::NMatrix> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Array2D_NMatrix__exposer.def( "__setstate__", &__setstate__base64< ::SireBase::Array2D<SireMaths::NMatrix> > );
        Array2D_NMatrix__exposer.def( "__getstate__", &__getstate__base64< ::SireBase::Array2D<SireMaths::NMatrix> > );
        Array2D_NMatrix__exposer.def( "__str__", &__str__< ::SireBase::Array2D<SireMaths::NMatrix> > );
        Array2D_NMatrix__exposer.def( "__repr__", &__str__< ::SireBase::Array2D<SireMaths::NMatrix> > );
    }

}
