// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "AngleRestraint.pypp.hpp"

namespace bp = boost::python;

#include "SireCAS/conditional.h"

#include "SireCAS/errors.h"

#include "SireCAS/power.h"

#include "SireCAS/symbols.h"

#include "SireCAS/values.h"

#include "SireFF/forcetable.h"

#include "SireID/index.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireUnits/angle.h"

#include "SireUnits/units.h"

#include "anglerestraint.h"

#include "anglerestraint.h"

SireMM::AngleRestraint __copy__(const SireMM::AngleRestraint &other){ return SireMM::AngleRestraint(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_AngleRestraint_class(){

    { //::SireMM::AngleRestraint
        typedef bp::class_< SireMM::AngleRestraint, bp::bases< SireMM::Restraint3D, SireMM::Restraint, SireBase::Property > > AngleRestraint_exposer_t;
        AngleRestraint_exposer_t AngleRestraint_exposer = AngleRestraint_exposer_t( "AngleRestraint", "This is a restraint that operates on the angle between\nthree SireMM::Point objects (e.g. three atoms in a molecule)\n\nAuthor: Christopher Woods\n", bp::init< >("Constructor") );
        bp::scope AngleRestraint_scope( AngleRestraint_exposer );
        AngleRestraint_exposer.def( bp::init< SireFF::PointRef const &, SireFF::PointRef const &, SireFF::PointRef const &, SireCAS::Expression const & >(( bp::arg("point0"), bp::arg("point1"), bp::arg("point2"), bp::arg("restraint") ), "Construct a restraint that acts on the angle within the\nthree points point0, point1 and point2 (theta == a(012)),\nrestraining the angle within these points using the expression\nrestraint") );
        AngleRestraint_exposer.def( bp::init< SireFF::PointRef const &, SireFF::PointRef const &, SireFF::PointRef const &, SireCAS::Expression const &, SireCAS::Values const & >(( bp::arg("point0"), bp::arg("point1"), bp::arg("point2"), bp::arg("restraint"), bp::arg("values") ), "Construct a restraint that acts on the angle within the\nthree points point0, point1 and point2 (theta == a(012)),\nrestraining the angle within these points using the expression\nrestraint") );
        AngleRestraint_exposer.def( bp::init< SireMM::AngleRestraint const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMM::AngleRestraint::builtinSymbols
        
            typedef ::SireCAS::Symbols ( ::SireMM::AngleRestraint::*builtinSymbols_function_type)(  ) const;
            builtinSymbols_function_type builtinSymbols_function_value( &::SireMM::AngleRestraint::builtinSymbols );
            
            AngleRestraint_exposer.def( 
                "builtinSymbols"
                , builtinSymbols_function_value
                , bp::release_gil_policy()
                , "Return the built-in symbols of this restraint" );
        
        }
        { //::SireMM::AngleRestraint::builtinValues
        
            typedef ::SireCAS::Values ( ::SireMM::AngleRestraint::*builtinValues_function_type)(  ) const;
            builtinValues_function_type builtinValues_function_value( &::SireMM::AngleRestraint::builtinValues );
            
            AngleRestraint_exposer.def( 
                "builtinValues"
                , builtinValues_function_value
                , bp::release_gil_policy()
                , "Return the values of the built-in symbols of this restraint" );
        
        }
        { //::SireMM::AngleRestraint::contains
        
            typedef bool ( ::SireMM::AngleRestraint::*contains_function_type)( ::SireMol::MolNum ) const;
            contains_function_type contains_function_value( &::SireMM::AngleRestraint::contains );
            
            AngleRestraint_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("molnum") )
                , bp::release_gil_policy()
                , "Return whether or not this restraint affects the molecule\nwith number molnum" );
        
        }
        { //::SireMM::AngleRestraint::contains
        
            typedef bool ( ::SireMM::AngleRestraint::*contains_function_type)( ::SireMol::MolID const & ) const;
            contains_function_type contains_function_value( &::SireMM::AngleRestraint::contains );
            
            AngleRestraint_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("molid") )
                , bp::release_gil_policy()
                , "Return whether or not this restraint affects the molecule\nwith ID molid" );
        
        }
        { //::SireMM::AngleRestraint::differentialRestraintFunction
        
            typedef ::SireCAS::Expression const & ( ::SireMM::AngleRestraint::*differentialRestraintFunction_function_type)(  ) const;
            differentialRestraintFunction_function_type differentialRestraintFunction_function_value( &::SireMM::AngleRestraint::differentialRestraintFunction );
            
            AngleRestraint_exposer.def( 
                "differentialRestraintFunction"
                , differentialRestraintFunction_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the function used to calculate the restraint force" );
        
        }
        { //::SireMM::AngleRestraint::differentiate
        
            typedef ::SireMM::RestraintPtr ( ::SireMM::AngleRestraint::*differentiate_function_type)( ::SireCAS::Symbol const & ) const;
            differentiate_function_type differentiate_function_value( &::SireMM::AngleRestraint::differentiate );
            
            AngleRestraint_exposer.def( 
                "differentiate"
                , differentiate_function_value
                , ( bp::arg("symbol") )
                , bp::release_gil_policy()
                , "Return the differential of this restraint with respect to\nthe symbol symbol\nThrow: SireCAS::unavailable_differential\n" );
        
        }
        { //::SireMM::AngleRestraint::force
        
            typedef void ( ::SireMM::AngleRestraint::*force_function_type)( ::SireFF::MolForceTable &,double ) const;
            force_function_type force_function_value( &::SireMM::AngleRestraint::force );
            
            AngleRestraint_exposer.def( 
                "force"
                , force_function_value
                , ( bp::arg("forcetable"), bp::arg("scale_force")=1 )
                , "Calculate the force acting on the molecule in the forcetable forcetable\ncaused by this restraint, and add it on to the forcetable scaled by\nscale_force" );
        
        }
        { //::SireMM::AngleRestraint::force
        
            typedef void ( ::SireMM::AngleRestraint::*force_function_type)( ::SireFF::ForceTable &,double ) const;
            force_function_type force_function_value( &::SireMM::AngleRestraint::force );
            
            AngleRestraint_exposer.def( 
                "force"
                , force_function_value
                , ( bp::arg("forcetable"), bp::arg("scale_force")=1 )
                , "Calculate the force acting on the molecules in the forcetable forcetable\ncaused by this restraint, and add it on to the forcetable scaled by\nscale_force" );
        
        }
        { //::SireMM::AngleRestraint::halfHarmonic
        
            typedef ::SireMM::AngleRestraint ( *halfHarmonic_function_type )( ::SireFF::PointRef const &,::SireFF::PointRef const &,::SireFF::PointRef const &,::SireUnits::Dimension::Angle const &,::SireMM::HarmonicAngleForceConstant const & );
            halfHarmonic_function_type halfHarmonic_function_value( &::SireMM::AngleRestraint::halfHarmonic );
            
            AngleRestraint_exposer.def( 
                "halfHarmonic"
                , halfHarmonic_function_value
                , ( bp::arg("point0"), bp::arg("point1"), bp::arg("point2"), bp::arg("angle"), bp::arg("force_constant") )
                , bp::release_gil_policy()
                , "Return a distance restraint that applied a half-harmonic potential\nbetween the points point0 and point1 above a distance distance\nusing a force constant force_constant" );
        
        }
        { //::SireMM::AngleRestraint::harmonic
        
            typedef ::SireMM::AngleRestraint ( *harmonic_function_type )( ::SireFF::PointRef const &,::SireFF::PointRef const &,::SireFF::PointRef const &,::SireMM::HarmonicAngleForceConstant const & );
            harmonic_function_type harmonic_function_value( &::SireMM::AngleRestraint::harmonic );
            
            AngleRestraint_exposer.def( 
                "harmonic"
                , harmonic_function_value
                , ( bp::arg("point0"), bp::arg("point1"), bp::arg("point2"), bp::arg("force_constant") )
                , bp::release_gil_policy()
                , "Return a distance restraint that applies a harmonic potential between\nthe points point0 and point1 using a force constant force_constant" );
        
        }
        { //::SireMM::AngleRestraint::molecules
        
            typedef ::SireMol::Molecules ( ::SireMM::AngleRestraint::*molecules_function_type)(  ) const;
            molecules_function_type molecules_function_value( &::SireMM::AngleRestraint::molecules );
            
            AngleRestraint_exposer.def( 
                "molecules"
                , molecules_function_value
                , bp::release_gil_policy()
                , "Return the molecules used in this restraint" );
        
        }
        { //::SireMM::AngleRestraint::nPoints
        
            typedef int ( ::SireMM::AngleRestraint::*nPoints_function_type)(  ) const;
            nPoints_function_type nPoints_function_value( &::SireMM::AngleRestraint::nPoints );
            
            AngleRestraint_exposer.def( 
                "nPoints"
                , nPoints_function_value
                , bp::release_gil_policy()
                , "This restraint involves three points" );
        
        }
        AngleRestraint_exposer.def( bp::self != bp::self );
        { //::SireMM::AngleRestraint::operator=
        
            typedef ::SireMM::AngleRestraint & ( ::SireMM::AngleRestraint::*assign_function_type)( ::SireMM::AngleRestraint const & ) ;
            assign_function_type assign_function_value( &::SireMM::AngleRestraint::operator= );
            
            AngleRestraint_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        AngleRestraint_exposer.def( bp::self == bp::self );
        { //::SireMM::AngleRestraint::point
        
            typedef ::SireFF::Point const & ( ::SireMM::AngleRestraint::*point_function_type)( int ) const;
            point_function_type point_function_value( &::SireMM::AngleRestraint::point );
            
            AngleRestraint_exposer.def( 
                "point"
                , point_function_value
                , ( bp::arg("i") )
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the ith point" );
        
        }
        { //::SireMM::AngleRestraint::point0
        
            typedef ::SireFF::Point const & ( ::SireMM::AngleRestraint::*point0_function_type)(  ) const;
            point0_function_type point0_function_value( &::SireMM::AngleRestraint::point0 );
            
            AngleRestraint_exposer.def( 
                "point0"
                , point0_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the first point" );
        
        }
        { //::SireMM::AngleRestraint::point1
        
            typedef ::SireFF::Point const & ( ::SireMM::AngleRestraint::*point1_function_type)(  ) const;
            point1_function_type point1_function_value( &::SireMM::AngleRestraint::point1 );
            
            AngleRestraint_exposer.def( 
                "point1"
                , point1_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the second point" );
        
        }
        { //::SireMM::AngleRestraint::point2
        
            typedef ::SireFF::Point const & ( ::SireMM::AngleRestraint::*point2_function_type)(  ) const;
            point2_function_type point2_function_value( &::SireMM::AngleRestraint::point2 );
            
            AngleRestraint_exposer.def( 
                "point2"
                , point2_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the third point" );
        
        }
        { //::SireMM::AngleRestraint::setSpace
        
            typedef void ( ::SireMM::AngleRestraint::*setSpace_function_type)( ::SireVol::Space const & ) ;
            setSpace_function_type setSpace_function_value( &::SireMM::AngleRestraint::setSpace );
            
            AngleRestraint_exposer.def( 
                "setSpace"
                , setSpace_function_value
                , ( bp::arg("space") )
                , bp::release_gil_policy()
                , "Set the space used to evaluate the energy of this restraint\nThrow: SireVol::incompatible_space\n" );
        
        }
        { //::SireMM::AngleRestraint::theta
        
            typedef ::SireCAS::Symbol const & ( *theta_function_type )(  );
            theta_function_type theta_function_value( &::SireMM::AngleRestraint::theta );
            
            AngleRestraint_exposer.def( 
                "theta"
                , theta_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the symbol that represents the angle between the points (theta)" );
        
        }
        { //::SireMM::AngleRestraint::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMM::AngleRestraint::typeName );
            
            AngleRestraint_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMM::AngleRestraint::update
        
            typedef void ( ::SireMM::AngleRestraint::*update_function_type)( ::SireMol::MoleculeData const & ) ;
            update_function_type update_function_value( &::SireMM::AngleRestraint::update );
            
            AngleRestraint_exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("moldata") )
                , bp::release_gil_policy()
                , "Update the points of this restraint using new molecule data from moldata\nThrow: SireBase::missing_property\nThrow: SireError::invalid_cast\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireMM::AngleRestraint::update
        
            typedef void ( ::SireMM::AngleRestraint::*update_function_type)( ::SireMol::Molecules const & ) ;
            update_function_type update_function_value( &::SireMM::AngleRestraint::update );
            
            AngleRestraint_exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("molecules") )
                , bp::release_gil_policy()
                , "Update the points of this restraint using new molecule data from molecules\nThrow: SireBase::missing_property\nThrow: SireError::invalid_cast\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireMM::AngleRestraint::usesMoleculesIn
        
            typedef bool ( ::SireMM::AngleRestraint::*usesMoleculesIn_function_type)( ::SireFF::ForceTable const & ) const;
            usesMoleculesIn_function_type usesMoleculesIn_function_value( &::SireMM::AngleRestraint::usesMoleculesIn );
            
            AngleRestraint_exposer.def( 
                "usesMoleculesIn"
                , usesMoleculesIn_function_value
                , ( bp::arg("forcetable") )
                , bp::release_gil_policy()
                , "Return whether or not this restraint involves any of the molecules\nthat are in the forcetable forcetable" );
        
        }
        { //::SireMM::AngleRestraint::usesMoleculesIn
        
            typedef bool ( ::SireMM::AngleRestraint::*usesMoleculesIn_function_type)( ::SireMol::Molecules const & ) const;
            usesMoleculesIn_function_type usesMoleculesIn_function_value( &::SireMM::AngleRestraint::usesMoleculesIn );
            
            AngleRestraint_exposer.def( 
                "usesMoleculesIn"
                , usesMoleculesIn_function_value
                , ( bp::arg("molecules") )
                , bp::release_gil_policy()
                , "Return whether or not this restraint involves any of the molecules\nin molecules" );
        
        }
        AngleRestraint_exposer.staticmethod( "halfHarmonic" );
        AngleRestraint_exposer.staticmethod( "harmonic" );
        AngleRestraint_exposer.staticmethod( "theta" );
        AngleRestraint_exposer.staticmethod( "typeName" );
        AngleRestraint_exposer.def( "__copy__", &__copy__);
        AngleRestraint_exposer.def( "__deepcopy__", &__copy__);
        AngleRestraint_exposer.def( "clone", &__copy__);
        AngleRestraint_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMM::AngleRestraint >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AngleRestraint_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMM::AngleRestraint >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AngleRestraint_exposer.def_pickle(sire_pickle_suite< ::SireMM::AngleRestraint >());
        AngleRestraint_exposer.def( "__str__", &__str__< ::SireMM::AngleRestraint > );
        AngleRestraint_exposer.def( "__repr__", &__str__< ::SireMM::AngleRestraint > );
    }

}
