// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "QMForce.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "qmmm.h"

#include "SireError/errors.h"

#include "qmmm.h"

const char* pvt_get_name(const SireOpenMM::QMForce&){ return "SireOpenMM::QMForce";}

#include "Helpers/release_gil_policy.hpp"

#include "SireError/errors.h"

#include "qmmm.h"

#include "SireError/errors.h"

#include "qmmm.h"

const char* pvt_get_name(const SireOpenMM::QMForce&){ return "SireOpenMM::QMForce";}

#include "Helpers/release_gil_policy.hpp"

void register_QMForce_class(){

    { //::SireOpenMM::QMForce
        typedef bp::class_< SireOpenMM::QMForce, boost::noncopyable > QMForce_exposer_t;
        QMForce_exposer_t QMForce_exposer = QMForce_exposer_t( "QMForce", "", bp::no_init );
        bp::scope QMForce_scope( QMForce_exposer );
        { //::SireOpenMM::QMForce::setLambda
        
            typedef void ( ::SireOpenMM::QMForce::*setLambda_function_type)( double ) ;
            setLambda_function_type setLambda_function_value( &::SireOpenMM::QMForce::setLambda );
            
            QMForce_exposer.def( 
                "setLambda"
                , setLambda_function_value
                , ( bp::arg("lambda") )
                , bp::release_gil_policy()
                , "Set the lambda weighting factor." );
        
        }
        QMForce_exposer.def( "__str__", &pvt_get_name);
        QMForce_exposer.def( "__repr__", &pvt_get_name);
        QMForce_exposer.def( "__str__", &pvt_get_name);
        QMForce_exposer.def( "__repr__", &pvt_get_name);
    }

}
