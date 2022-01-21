// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "AngleComponent.pypp.hpp"

namespace bp = boost::python;

#include "SireID/index.h"

#include "SireMaths/triangle.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireSystem/system.h"

#include "SireUnits/units.h"

#include "SireVol/space.h"

#include "anglecomponent.h"

#include "delta.h"

#include "anglecomponent.h"

SireSystem::AngleComponent __copy__(const SireSystem::AngleComponent &other){ return SireSystem::AngleComponent(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_AngleComponent_class(){

    { //::SireSystem::AngleComponent
        typedef bp::class_< SireSystem::AngleComponent, bp::bases< SireSystem::GeometryComponent, SireSystem::Constraint, SireBase::Property > > AngleComponent_exposer_t;
        AngleComponent_exposer_t AngleComponent_exposer = AngleComponent_exposer_t( "AngleComponent", "This is a constraint that constrains a symbol to equal the\nvalue of an expression that involves an angle between three points\nor atoms\n\nAuthor: Christopher Woods\n", bp::init< >("Null constructor") );
        bp::scope AngleComponent_scope( AngleComponent_exposer );
        AngleComponent_exposer.def( bp::init< SireCAS::Symbol const &, SireFF::PointRef const &, SireFF::PointRef const &, SireFF::PointRef const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("constrained_symbol"), bp::arg("point0"), bp::arg("point1"), bp::arg("point2"), bp::arg("map")=SireBase::PropertyMap() ), "Construct to set the value of constrained_symbol equal to the\nangle between the three points point0, point1 and point2") );
        AngleComponent_exposer.def( bp::init< SireCAS::Symbol const &, SireFF::PointRef const &, SireFF::PointRef const &, SireFF::PointRef const &, SireCAS::Expression const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("constrained_symbol"), bp::arg("point0"), bp::arg("point1"), bp::arg("point2"), bp::arg("geometry_expression"), bp::arg("map")=SireBase::PropertyMap() ), "Construct to set the value of constrained_symbol equal to the\nexpression based on the angles within the three points\npoint0, point1 and point2") );
        AngleComponent_exposer.def( bp::init< SireSystem::AngleComponent const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireSystem::AngleComponent::nPoints
        
            typedef int ( ::SireSystem::AngleComponent::*nPoints_function_type)(  ) const;
            nPoints_function_type nPoints_function_value( &::SireSystem::AngleComponent::nPoints );
            
            AngleComponent_exposer.def( 
                "nPoints"
                , nPoints_function_value
                , "Return the number of points (3)" );
        
        }
        AngleComponent_exposer.def( bp::self != bp::self );
        { //::SireSystem::AngleComponent::operator=
        
            typedef ::SireSystem::AngleComponent & ( ::SireSystem::AngleComponent::*assign_function_type)( ::SireSystem::AngleComponent const & ) ;
            assign_function_type assign_function_value( &::SireSystem::AngleComponent::operator= );
            
            AngleComponent_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        AngleComponent_exposer.def( bp::self == bp::self );
        { //::SireSystem::AngleComponent::point
        
            typedef ::SireFF::Point const & ( ::SireSystem::AngleComponent::*point_function_type)( int ) const;
            point_function_type point_function_value( &::SireSystem::AngleComponent::point );
            
            AngleComponent_exposer.def( 
                "point"
                , point_function_value
                , ( bp::arg("i") )
                , bp::return_value_policy<bp::clone_const_reference>()
                , "Return the ith point\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireSystem::AngleComponent::point0
        
            typedef ::SireFF::Point const & ( ::SireSystem::AngleComponent::*point0_function_type)(  ) const;
            point0_function_type point0_function_value( &::SireSystem::AngleComponent::point0 );
            
            AngleComponent_exposer.def( 
                "point0"
                , point0_function_value
                , bp::return_value_policy<bp::clone_const_reference>()
                , "Return the first point between which the angle is calculated" );
        
        }
        { //::SireSystem::AngleComponent::point1
        
            typedef ::SireFF::Point const & ( ::SireSystem::AngleComponent::*point1_function_type)(  ) const;
            point1_function_type point1_function_value( &::SireSystem::AngleComponent::point1 );
            
            AngleComponent_exposer.def( 
                "point1"
                , point1_function_value
                , bp::return_value_policy<bp::clone_const_reference>()
                , "Return the second point between which the angle is calculated" );
        
        }
        { //::SireSystem::AngleComponent::point2
        
            typedef ::SireFF::Point const & ( ::SireSystem::AngleComponent::*point2_function_type)(  ) const;
            point2_function_type point2_function_value( &::SireSystem::AngleComponent::point2 );
            
            AngleComponent_exposer.def( 
                "point2"
                , point2_function_value
                , bp::return_value_policy<bp::clone_const_reference>()
                , "Return the third point between which the angle is calculated" );
        
        }
        { //::SireSystem::AngleComponent::r01
        
            typedef ::SireCAS::Symbol const & ( *r01_function_type )(  );
            r01_function_type r01_function_value( &::SireSystem::AngleComponent::r01 );
            
            AngleComponent_exposer.def( 
                "r01"
                , r01_function_value
                , bp::return_value_policy<bp::clone_const_reference>()
                , "Return the symbol that represents the 0-1 distance" );
        
        }
        { //::SireSystem::AngleComponent::r02
        
            typedef ::SireCAS::Symbol const & ( *r02_function_type )(  );
            r02_function_type r02_function_value( &::SireSystem::AngleComponent::r02 );
            
            AngleComponent_exposer.def( 
                "r02"
                , r02_function_value
                , bp::return_value_policy<bp::clone_const_reference>()
                , "Return the symbol that represents the 0-2 distance" );
        
        }
        { //::SireSystem::AngleComponent::r12
        
            typedef ::SireCAS::Symbol const & ( *r12_function_type )(  );
            r12_function_type r12_function_value( &::SireSystem::AngleComponent::r12 );
            
            AngleComponent_exposer.def( 
                "r12"
                , r12_function_value
                , bp::return_value_policy<bp::clone_const_reference>()
                , "Return the symbol that represents the 1-2 distance" );
        
        }
        { //::SireSystem::AngleComponent::theta
        
            typedef ::SireCAS::Symbol const & ( *theta_function_type )(  );
            theta_function_type theta_function_value( &::SireSystem::AngleComponent::theta );
            
            AngleComponent_exposer.def( 
                "theta"
                , theta_function_value
                , bp::return_value_policy<bp::clone_const_reference>()
                , "Return the symbol that represents the central angle between\nthe three points" );
        
        }
        { //::SireSystem::AngleComponent::theta012
        
            typedef ::SireCAS::Symbol const & ( *theta012_function_type )(  );
            theta012_function_type theta012_function_value( &::SireSystem::AngleComponent::theta012 );
            
            AngleComponent_exposer.def( 
                "theta012"
                , theta012_function_value
                , bp::return_value_policy<bp::clone_const_reference>()
                , "Return the symbol that represents the angle 012 between\nthe three points" );
        
        }
        { //::SireSystem::AngleComponent::theta021
        
            typedef ::SireCAS::Symbol const & ( *theta021_function_type )(  );
            theta021_function_type theta021_function_value( &::SireSystem::AngleComponent::theta021 );
            
            AngleComponent_exposer.def( 
                "theta021"
                , theta021_function_value
                , bp::return_value_policy<bp::clone_const_reference>()
                , "Return the symbol that represents the angle 021 between\nthe three points" );
        
        }
        { //::SireSystem::AngleComponent::theta102
        
            typedef ::SireCAS::Symbol const & ( *theta102_function_type )(  );
            theta102_function_type theta102_function_value( &::SireSystem::AngleComponent::theta102 );
            
            AngleComponent_exposer.def( 
                "theta102"
                , theta102_function_value
                , bp::return_value_policy<bp::clone_const_reference>()
                , "Return the symbol that represents the angle 102 between\nthe three points" );
        
        }
        { //::SireSystem::AngleComponent::toString
        
            typedef ::QString ( ::SireSystem::AngleComponent::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireSystem::AngleComponent::toString );
            
            AngleComponent_exposer.def( 
                "toString"
                , toString_function_value
                , "" );
        
        }
        { //::SireSystem::AngleComponent::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireSystem::AngleComponent::typeName );
            
            AngleComponent_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        AngleComponent_exposer.staticmethod( "r01" );
        AngleComponent_exposer.staticmethod( "r02" );
        AngleComponent_exposer.staticmethod( "r12" );
        AngleComponent_exposer.staticmethod( "theta" );
        AngleComponent_exposer.staticmethod( "theta012" );
        AngleComponent_exposer.staticmethod( "theta021" );
        AngleComponent_exposer.staticmethod( "theta102" );
        AngleComponent_exposer.staticmethod( "typeName" );
        AngleComponent_exposer.def( "__copy__", &__copy__);
        AngleComponent_exposer.def( "__deepcopy__", &__copy__);
        AngleComponent_exposer.def( "clone", &__copy__);
        AngleComponent_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireSystem::AngleComponent >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AngleComponent_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireSystem::AngleComponent >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AngleComponent_exposer.def( "__setstate__", &__setstate__base64< ::SireSystem::AngleComponent > );
        AngleComponent_exposer.def( "__getstate__", &__getstate__base64< ::SireSystem::AngleComponent > );
        AngleComponent_exposer.def( "__str__", &__str__< ::SireSystem::AngleComponent > );
        AngleComponent_exposer.def( "__repr__", &__str__< ::SireSystem::AngleComponent > );
    }

}
