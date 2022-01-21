// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "SwitchingFunction.pypp.hpp"

namespace bp = boost::python;

#include "SireFF/errors.h"

#include "SireMaths/maths.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireUnits/units.h"

#include "switchingfunction.h"

#include <QMutex>

#include <cmath>

#include <numeric>

#include "switchingfunction.h"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_SwitchingFunction_class(){

    { //::SireMM::SwitchingFunction
        typedef bp::class_< SireMM::SwitchingFunction, bp::bases< SireBase::Property >, boost::noncopyable > SwitchingFunction_exposer_t;
        SwitchingFunction_exposer_t SwitchingFunction_exposer = SwitchingFunction_exposer_t( "SwitchingFunction", "This is the virtual base class of all switching functions. These\nreturn scale factors based on the supplied distance\n\nAuthor: Christopher Woods\n", bp::no_init );
        bp::scope SwitchingFunction_scope( SwitchingFunction_exposer );
        { //::SireMM::SwitchingFunction::cutoffDistance
        
            typedef ::SireUnits::Dimension::Length ( ::SireMM::SwitchingFunction::*cutoffDistance_function_type)(  ) const;
            cutoffDistance_function_type cutoffDistance_function_value( &::SireMM::SwitchingFunction::cutoffDistance );
            
            SwitchingFunction_exposer.def( 
                "cutoffDistance"
                , cutoffDistance_function_value
                , "" );
        
        }
        { //::SireMM::SwitchingFunction::dElectrostaticScaleFactor
        
            typedef double ( ::SireMM::SwitchingFunction::*dElectrostaticScaleFactor_function_type)( ::SireUnits::Dimension::Length ) const;
            dElectrostaticScaleFactor_function_type dElectrostaticScaleFactor_function_value( &::SireMM::SwitchingFunction::dElectrostaticScaleFactor );
            
            SwitchingFunction_exposer.def( 
                "dElectrostaticScaleFactor"
                , dElectrostaticScaleFactor_function_value
                , ( bp::arg("dist") )
                , "Return the derivative (gradient) of the electrostatic\nscale factor at the distance dist" );
        
        }
        { //::SireMM::SwitchingFunction::dVDWScaleFactor
        
            typedef double ( ::SireMM::SwitchingFunction::*dVDWScaleFactor_function_type)( ::SireUnits::Dimension::Length ) const;
            dVDWScaleFactor_function_type dVDWScaleFactor_function_value( &::SireMM::SwitchingFunction::dVDWScaleFactor );
            
            SwitchingFunction_exposer.def( 
                "dVDWScaleFactor"
                , dVDWScaleFactor_function_value
                , ( bp::arg("dist") )
                , "Return the derivative (gradient) of the VDW\nscale factor at the distance dist" );
        
        }
        { //::SireMM::SwitchingFunction::electrostaticCutoffDistance
        
            typedef ::SireUnits::Dimension::Length ( ::SireMM::SwitchingFunction::*electrostaticCutoffDistance_function_type)(  ) const;
            electrostaticCutoffDistance_function_type electrostaticCutoffDistance_function_value( &::SireMM::SwitchingFunction::electrostaticCutoffDistance );
            
            SwitchingFunction_exposer.def( 
                "electrostaticCutoffDistance"
                , electrostaticCutoffDistance_function_value
                , "" );
        
        }
        { //::SireMM::SwitchingFunction::electrostaticFeatherDistance
        
            typedef ::SireUnits::Dimension::Length ( ::SireMM::SwitchingFunction::*electrostaticFeatherDistance_function_type)(  ) const;
            electrostaticFeatherDistance_function_type electrostaticFeatherDistance_function_value( &::SireMM::SwitchingFunction::electrostaticFeatherDistance );
            
            SwitchingFunction_exposer.def( 
                "electrostaticFeatherDistance"
                , electrostaticFeatherDistance_function_value
                , "" );
        
        }
        { //::SireMM::SwitchingFunction::electrostaticScaleFactor
        
            typedef double ( ::SireMM::SwitchingFunction::*electrostaticScaleFactor_function_type)( ::SireUnits::Dimension::Length ) const;
            electrostaticScaleFactor_function_type electrostaticScaleFactor_function_value( &::SireMM::SwitchingFunction::electrostaticScaleFactor );
            
            SwitchingFunction_exposer.def( 
                "electrostaticScaleFactor"
                , electrostaticScaleFactor_function_value
                , ( bp::arg("dist") )
                , "Return the electrostatic scale factor for the distance dist" );
        
        }
        { //::SireMM::SwitchingFunction::featherDistance
        
            typedef ::SireUnits::Dimension::Length ( ::SireMM::SwitchingFunction::*featherDistance_function_type)(  ) const;
            featherDistance_function_type featherDistance_function_value( &::SireMM::SwitchingFunction::featherDistance );
            
            SwitchingFunction_exposer.def( 
                "featherDistance"
                , featherDistance_function_value
                , "" );
        
        }
        { //::SireMM::SwitchingFunction::null
        
            typedef ::SireMM::NoCutoff const & ( *null_function_type )(  );
            null_function_type null_function_value( &::SireMM::SwitchingFunction::null );
            
            SwitchingFunction_exposer.def( 
                "null"
                , null_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMM::SwitchingFunction::toString
        
            typedef ::QString ( ::SireMM::SwitchingFunction::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMM::SwitchingFunction::toString );
            
            SwitchingFunction_exposer.def( 
                "toString"
                , toString_function_value
                , "Return a string representation of this switching function" );
        
        }
        { //::SireMM::SwitchingFunction::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMM::SwitchingFunction::typeName );
            
            SwitchingFunction_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireMM::SwitchingFunction::vdwCutoffDistance
        
            typedef ::SireUnits::Dimension::Length ( ::SireMM::SwitchingFunction::*vdwCutoffDistance_function_type)(  ) const;
            vdwCutoffDistance_function_type vdwCutoffDistance_function_value( &::SireMM::SwitchingFunction::vdwCutoffDistance );
            
            SwitchingFunction_exposer.def( 
                "vdwCutoffDistance"
                , vdwCutoffDistance_function_value
                , "" );
        
        }
        { //::SireMM::SwitchingFunction::vdwFeatherDistance
        
            typedef ::SireUnits::Dimension::Length ( ::SireMM::SwitchingFunction::*vdwFeatherDistance_function_type)(  ) const;
            vdwFeatherDistance_function_type vdwFeatherDistance_function_value( &::SireMM::SwitchingFunction::vdwFeatherDistance );
            
            SwitchingFunction_exposer.def( 
                "vdwFeatherDistance"
                , vdwFeatherDistance_function_value
                , "" );
        
        }
        { //::SireMM::SwitchingFunction::vdwScaleFactor
        
            typedef double ( ::SireMM::SwitchingFunction::*vdwScaleFactor_function_type)( ::SireUnits::Dimension::Length ) const;
            vdwScaleFactor_function_type vdwScaleFactor_function_value( &::SireMM::SwitchingFunction::vdwScaleFactor );
            
            SwitchingFunction_exposer.def( 
                "vdwScaleFactor"
                , vdwScaleFactor_function_value
                , ( bp::arg("dist") )
                , "Return the VDW scale factor for the distance dist" );
        
        }
        SwitchingFunction_exposer.staticmethod( "null" );
        SwitchingFunction_exposer.staticmethod( "typeName" );
        SwitchingFunction_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMM::SwitchingFunction >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        SwitchingFunction_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMM::SwitchingFunction >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        SwitchingFunction_exposer.def( "__setstate__", &__setstate__base64< ::SireMM::SwitchingFunction > );
        SwitchingFunction_exposer.def( "__getstate__", &__getstate__base64< ::SireMM::SwitchingFunction > );
        SwitchingFunction_exposer.def( "__str__", &__str__< ::SireMM::SwitchingFunction > );
        SwitchingFunction_exposer.def( "__repr__", &__str__< ::SireMM::SwitchingFunction > );
    }

}
