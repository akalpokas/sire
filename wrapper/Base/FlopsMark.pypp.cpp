// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "FlopsMark.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "countflops.h"

#include <QDebug>

#include <cmath>

#include "countflops.h"

SireBase::FlopsMark __copy__(const SireBase::FlopsMark &other){ return SireBase::FlopsMark(other); }

const char* pvt_get_name(const SireBase::FlopsMark&){ return "SireBase::FlopsMark";}

void register_FlopsMark_class(){

    { //::SireBase::FlopsMark
        typedef bp::class_< SireBase::FlopsMark > FlopsMark_exposer_t;
        FlopsMark_exposer_t FlopsMark_exposer = FlopsMark_exposer_t( "FlopsMark", bp::init< >() );
        bp::scope FlopsMark_scope( FlopsMark_exposer );
        FlopsMark_exposer.def( bp::init< SireBase::FlopsMark const & >(( bp::arg("other") )) );
        { //::SireBase::FlopsMark::benchmark
        
            typedef double ( *benchmark_function_type )(  );
            benchmark_function_type benchmark_function_value( &::SireBase::FlopsMark::benchmark );
            
            FlopsMark_exposer.def( 
                "benchmark"
                , benchmark_function_value );
        
        }
        { //::SireBase::FlopsMark::benchmarkProduct
        
            typedef double ( *benchmarkProduct_function_type )(  );
            benchmarkProduct_function_type benchmarkProduct_function_value( &::SireBase::FlopsMark::benchmarkProduct );
            
            FlopsMark_exposer.def( 
                "benchmarkProduct"
                , benchmarkProduct_function_value );
        
        }
        { //::SireBase::FlopsMark::benchmarkQuotient
        
            typedef double ( *benchmarkQuotient_function_type )(  );
            benchmarkQuotient_function_type benchmarkQuotient_function_value( &::SireBase::FlopsMark::benchmarkQuotient );
            
            FlopsMark_exposer.def( 
                "benchmarkQuotient"
                , benchmarkQuotient_function_value );
        
        }
        { //::SireBase::FlopsMark::benchmarkSum
        
            typedef double ( *benchmarkSum_function_type )(  );
            benchmarkSum_function_type benchmarkSum_function_value( &::SireBase::FlopsMark::benchmarkSum );
            
            FlopsMark_exposer.def( 
                "benchmarkSum"
                , benchmarkSum_function_value );
        
        }
        { //::SireBase::FlopsMark::nFlops
        
            typedef int ( ::SireBase::FlopsMark::*nFlops_function_type)(  ) const;
            nFlops_function_type nFlops_function_value( &::SireBase::FlopsMark::nFlops );
            
            FlopsMark_exposer.def( 
                "nFlops"
                , nFlops_function_value );
        
        }
        { //::SireBase::FlopsMark::nThreads
        
            typedef int ( ::SireBase::FlopsMark::*nThreads_function_type)(  ) const;
            nThreads_function_type nThreads_function_value( &::SireBase::FlopsMark::nThreads );
            
            FlopsMark_exposer.def( 
                "nThreads"
                , nThreads_function_value );
        
        }
        FlopsMark_exposer.def( bp::self - bp::self );
        { //::SireBase::FlopsMark::operator=
        
            typedef ::SireBase::FlopsMark & ( ::SireBase::FlopsMark::*assign_function_type)( ::SireBase::FlopsMark const & ) ;
            assign_function_type assign_function_value( &::SireBase::FlopsMark::operator= );
            
            FlopsMark_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        { //::SireBase::FlopsMark::operator[]
        
            typedef ::SireBase::FlopsMark ( ::SireBase::FlopsMark::*__getitem___function_type)( int ) const;
            __getitem___function_type __getitem___function_value( &::SireBase::FlopsMark::operator[] );
            
            FlopsMark_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") ) );
        
        }
        { //::SireBase::FlopsMark::threadFlops
        
            typedef ::SireBase::FlopsMark ( ::SireBase::FlopsMark::*threadFlops_function_type)( int ) const;
            threadFlops_function_type threadFlops_function_value( &::SireBase::FlopsMark::threadFlops );
            
            FlopsMark_exposer.def( 
                "threadFlops"
                , threadFlops_function_value
                , ( bp::arg("i") ) );
        
        }
        FlopsMark_exposer.staticmethod( "benchmark" );
        FlopsMark_exposer.staticmethod( "benchmarkProduct" );
        FlopsMark_exposer.staticmethod( "benchmarkQuotient" );
        FlopsMark_exposer.staticmethod( "benchmarkSum" );
        FlopsMark_exposer.def( "__copy__", &__copy__);
        FlopsMark_exposer.def( "__deepcopy__", &__copy__);
        FlopsMark_exposer.def( "clone", &__copy__);
        FlopsMark_exposer.def( "__str__", &pvt_get_name);
        FlopsMark_exposer.def( "__repr__", &pvt_get_name);
    }

}
