// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "SoftCLJPotentialInterface_InterSoftCLJPotential_.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/mover.hpp"

#include "SireMol/partialmolecule.h"

#include "intersoftcljff.h"

#include "intersoftcljff.h"

const char* pvt_get_name(const SireMM::SoftCLJPotentialInterface<SireMM::InterSoftCLJPotential>&){ return "SireMM::SoftCLJPotentialInterface<SireMM::InterSoftCLJPotential>";}

#include "Helpers/release_gil_policy.hpp"

void register_SoftCLJPotentialInterface_InterSoftCLJPotential__class(){

    { //::SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential >
        typedef bp::class_< SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential >, bp::bases< SireMM::CLJPotentialInterface<SireMM::InterSoftCLJPotential> >, boost::noncopyable > SoftCLJPotentialInterface_InterSoftCLJPotential__exposer_t;
        SoftCLJPotentialInterface_InterSoftCLJPotential__exposer_t SoftCLJPotentialInterface_InterSoftCLJPotential__exposer = SoftCLJPotentialInterface_InterSoftCLJPotential__exposer_t( "SoftCLJPotentialInterface_InterSoftCLJPotential_", "", bp::no_init );
        bp::scope SoftCLJPotentialInterface_InterSoftCLJPotential__scope( SoftCLJPotentialInterface_InterSoftCLJPotential__exposer );
        { //::SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential >::alpha
        
            typedef SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > exported_class_t;
            typedef double ( ::SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential >::*alpha_function_type)(  ) const;
            alpha_function_type alpha_function_value( &::SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential >::alpha );
            
            SoftCLJPotentialInterface_InterSoftCLJPotential__exposer.def( 
                "alpha"
                , alpha_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential >::alpha
        
            typedef SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > exported_class_t;
            typedef double ( ::SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential >::*alpha_function_type)( int ) const;
            alpha_function_type alpha_function_value( &::SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential >::alpha );
            
            SoftCLJPotentialInterface_InterSoftCLJPotential__exposer.def( 
                "alpha"
                , alpha_function_value
                , ( bp::arg("i") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential >::clearAlphas
        
            typedef SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > exported_class_t;
            typedef void ( ::SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential >::*clearAlphas_function_type)(  ) ;
            clearAlphas_function_type clearAlphas_function_value( &::SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential >::clearAlphas );
            
            SoftCLJPotentialInterface_InterSoftCLJPotential__exposer.def( 
                "clearAlphas"
                , clearAlphas_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential >::coulombPower
        
            typedef SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > exported_class_t;
            typedef int ( ::SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential >::*coulombPower_function_type)(  ) const;
            coulombPower_function_type coulombPower_function_value( &::SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential >::coulombPower );
            
            SoftCLJPotentialInterface_InterSoftCLJPotential__exposer.def( 
                "coulombPower"
                , coulombPower_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential >::hasAlphaValue
        
            typedef SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential >::*hasAlphaValue_function_type)( int ) const;
            hasAlphaValue_function_type hasAlphaValue_function_value( &::SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential >::hasAlphaValue );
            
            SoftCLJPotentialInterface_InterSoftCLJPotential__exposer.def( 
                "hasAlphaValue"
                , hasAlphaValue_function_value
                , ( bp::arg("i") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential >::ljPower
        
            typedef SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > exported_class_t;
            typedef int ( ::SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential >::*ljPower_function_type)(  ) const;
            ljPower_function_type ljPower_function_value( &::SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential >::ljPower );
            
            SoftCLJPotentialInterface_InterSoftCLJPotential__exposer.def( 
                "ljPower"
                , ljPower_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential >::nActiveAlphaComponents
        
            typedef SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > exported_class_t;
            typedef int ( ::SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential >::*nActiveAlphaComponents_function_type)(  ) const;
            nActiveAlphaComponents_function_type nActiveAlphaComponents_function_value( &::SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential >::nActiveAlphaComponents );
            
            SoftCLJPotentialInterface_InterSoftCLJPotential__exposer.def( 
                "nActiveAlphaComponents"
                , nActiveAlphaComponents_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential >::removeAlpha
        
            typedef SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential >::*removeAlpha_function_type)( int ) ;
            removeAlpha_function_type removeAlpha_function_value( &::SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential >::removeAlpha );
            
            SoftCLJPotentialInterface_InterSoftCLJPotential__exposer.def( 
                "removeAlpha"
                , removeAlpha_function_value
                , ( bp::arg("i") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential >::setAlpha
        
            typedef SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential >::*setAlpha_function_type)( double ) ;
            setAlpha_function_type setAlpha_function_value( &::SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential >::setAlpha );
            
            SoftCLJPotentialInterface_InterSoftCLJPotential__exposer.def( 
                "setAlpha"
                , setAlpha_function_value
                , ( bp::arg("alpha") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential >::setAlpha
        
            typedef SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential >::*setAlpha_function_type)( int,double ) ;
            setAlpha_function_type setAlpha_function_value( &::SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential >::setAlpha );
            
            SoftCLJPotentialInterface_InterSoftCLJPotential__exposer.def( 
                "setAlpha"
                , setAlpha_function_value
                , ( bp::arg("i"), bp::arg("alpha") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential >::setCoulombPower
        
            typedef SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential >::*setCoulombPower_function_type)( int ) ;
            setCoulombPower_function_type setCoulombPower_function_value( &::SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential >::setCoulombPower );
            
            SoftCLJPotentialInterface_InterSoftCLJPotential__exposer.def( 
                "setCoulombPower"
                , setCoulombPower_function_value
                , ( bp::arg("power") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential >::setLJPower
        
            typedef SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential >::*setLJPower_function_type)( int ) ;
            setLJPower_function_type setLJPower_function_value( &::SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential >::setLJPower );
            
            SoftCLJPotentialInterface_InterSoftCLJPotential__exposer.def( 
                "setLJPower"
                , setLJPower_function_value
                , ( bp::arg("power") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential >::setShiftDelta
        
            typedef SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > exported_class_t;
            typedef bool ( ::SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential >::*setShiftDelta_function_type)( double ) ;
            setShiftDelta_function_type setShiftDelta_function_value( &::SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential >::setShiftDelta );
            
            SoftCLJPotentialInterface_InterSoftCLJPotential__exposer.def( 
                "setShiftDelta"
                , setShiftDelta_function_value
                , ( bp::arg("delta") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential >::shiftDelta
        
            typedef SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential > exported_class_t;
            typedef double ( ::SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential >::*shiftDelta_function_type)(  ) const;
            shiftDelta_function_type shiftDelta_function_value( &::SireMM::SoftCLJPotentialInterface< SireMM::InterSoftCLJPotential >::shiftDelta );
            
            SoftCLJPotentialInterface_InterSoftCLJPotential__exposer.def( 
                "shiftDelta"
                , shiftDelta_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        SoftCLJPotentialInterface_InterSoftCLJPotential__exposer.def( "__str__", &pvt_get_name);
        SoftCLJPotentialInterface_InterSoftCLJPotential__exposer.def( "__repr__", &pvt_get_name);
    }

}
