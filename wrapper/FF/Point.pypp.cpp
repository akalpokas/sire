// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "Point.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/evaluator.h"

#include "SireMol/mgidx.h"

#include "SireMol/molecule.h"

#include "SireMol/moleculegroup.h"

#include "SireMol/moleculegroups.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireVol/aabox.h"

#include "SireVol/errors.h"

#include "forcetable.h"

#include "point.h"

#include <boost/tuple/tuple.hpp>

#include "point.h"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_Point_class(){

    { //::SireFF::Point
        typedef bp::class_< SireFF::Point, bp::bases< SireBase::Property >, boost::noncopyable > Point_exposer_t;
        Point_exposer_t Point_exposer = Point_exposer_t( "Point", "This is the base class of all Points. A Point is a class that\nallows a view of a molecule (or molecules) to be turned\ninto a 3D point. This 3D point can then be used within a restraint\nclass (thereby seperating out the code that implements the\nrestraint from the code that selects the parts of the molecule(s)\nbeing restrained)\n\nAuthor: Christopher Woods\n", bp::no_init );
        bp::scope Point_scope( Point_exposer );
        { //::SireFF::Point::addForce
        
            typedef bool ( ::SireFF::Point::*addForce_function_type)( ::SireFF::MolForceTable &,::SireMaths::Vector const & ) const;
            addForce_function_type addForce_function_value( &::SireFF::Point::addForce );
            
            Point_exposer.def( 
                "addForce"
                , addForce_function_value
                , ( bp::arg("molforces"), bp::arg("force") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Point::addForce
        
            typedef bool ( ::SireFF::Point::*addForce_function_type)( ::SireFF::ForceTable &,::SireMaths::Vector const & ) const;
            addForce_function_type addForce_function_value( &::SireFF::Point::addForce );
            
            Point_exposer.def( 
                "addForce"
                , addForce_function_value
                , ( bp::arg("forces"), bp::arg("force") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Point::areIntraMoleculePoints
        
            typedef bool ( *areIntraMoleculePoints_function_type )( ::SireFF::Point const &,::SireFF::Point const & );
            areIntraMoleculePoints_function_type areIntraMoleculePoints_function_value( &::SireFF::Point::areIntraMoleculePoints );
            
            Point_exposer.def( 
                "areIntraMoleculePoints"
                , areIntraMoleculePoints_function_value
                , ( bp::arg("point0"), bp::arg("point1") )
                , bp::release_gil_policy()
                , "Return whether or not the points point0 and point1 are\nboth within the same molecule (so together are intra-molecule points)" );
        
        }
        { //::SireFF::Point::areIntraMoleculePoints
        
            typedef bool ( *areIntraMoleculePoints_function_type )( ::SireFF::Point const &,::SireFF::Point const &,::SireFF::Point const & );
            areIntraMoleculePoints_function_type areIntraMoleculePoints_function_value( &::SireFF::Point::areIntraMoleculePoints );
            
            Point_exposer.def( 
                "areIntraMoleculePoints"
                , areIntraMoleculePoints_function_value
                , ( bp::arg("point0"), bp::arg("point1"), bp::arg("point2") )
                , bp::release_gil_policy()
                , "Return whether or not the points point0, point1 and point2 are\nwithin the same molecule (so together are intra-molecule points)" );
        
        }
        { //::SireFF::Point::areIntraMoleculePoints
        
            typedef bool ( *areIntraMoleculePoints_function_type )( ::SireFF::Point const &,::SireFF::Point const &,::SireFF::Point const &,::SireFF::Point const & );
            areIntraMoleculePoints_function_type areIntraMoleculePoints_function_value( &::SireFF::Point::areIntraMoleculePoints );
            
            Point_exposer.def( 
                "areIntraMoleculePoints"
                , areIntraMoleculePoints_function_value
                , ( bp::arg("point0"), bp::arg("point1"), bp::arg("point2"), bp::arg("point3") )
                , bp::release_gil_policy()
                , "Return whether or not the points point0, point1, point2 and point3 are\nwithin the same molecule (so together are intra-molecule points)" );
        
        }
        { //::SireFF::Point::areIntraMoleculePoints
        
            typedef bool ( *areIntraMoleculePoints_function_type )( ::QVector< SireBase::PropPtr< SireFF::Point > > const & );
            areIntraMoleculePoints_function_type areIntraMoleculePoints_function_value( &::SireFF::Point::areIntraMoleculePoints );
            
            Point_exposer.def( 
                "areIntraMoleculePoints"
                , areIntraMoleculePoints_function_value
                , ( bp::arg("points") )
                , bp::release_gil_policy()
                , "Return whether or not the points in points are all\nwithin the same molecule (so together are intra-molecule points)" );
        
        }
        { //::SireFF::Point::contains
        
            typedef bool ( ::SireFF::Point::*contains_function_type)( ::SireMol::MolNum ) const;
            contains_function_type contains_function_value( &::SireFF::Point::contains );
            
            Point_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("molnum") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Point::contains
        
            typedef bool ( ::SireFF::Point::*contains_function_type)( ::SireMol::MolID const & ) const;
            contains_function_type contains_function_value( &::SireFF::Point::contains );
            
            Point_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("molid") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Point::isExtraMoleculePoint
        
            typedef bool ( ::SireFF::Point::*isExtraMoleculePoint_function_type)(  ) const;
            isExtraMoleculePoint_function_type isExtraMoleculePoint_function_value( &::SireFF::Point::isExtraMoleculePoint );
            
            Point_exposer.def( 
                "isExtraMoleculePoint"
                , isExtraMoleculePoint_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Point::isInterMoleculePoint
        
            typedef bool ( ::SireFF::Point::*isInterMoleculePoint_function_type)(  ) const;
            isInterMoleculePoint_function_type isInterMoleculePoint_function_value( &::SireFF::Point::isInterMoleculePoint );
            
            Point_exposer.def( 
                "isInterMoleculePoint"
                , isInterMoleculePoint_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Point::isIntraMoleculePoint
        
            typedef bool ( ::SireFF::Point::*isIntraMoleculePoint_function_type)(  ) const;
            isIntraMoleculePoint_function_type isIntraMoleculePoint_function_value( &::SireFF::Point::isIntraMoleculePoint );
            
            Point_exposer.def( 
                "isIntraMoleculePoint"
                , isIntraMoleculePoint_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Point::molecules
        
            typedef ::SireMol::Molecules ( ::SireFF::Point::*molecules_function_type)(  ) const;
            molecules_function_type molecules_function_value( &::SireFF::Point::molecules );
            
            Point_exposer.def( 
                "molecules"
                , molecules_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Point::nMolecules
        
            typedef int ( ::SireFF::Point::*nMolecules_function_type)(  ) const;
            nMolecules_function_type nMolecules_function_value( &::SireFF::Point::nMolecules );
            
            Point_exposer.def( 
                "nMolecules"
                , nMolecules_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Point::null
        
            typedef ::SireFF::VectorPoint const & ( *null_function_type )(  );
            null_function_type null_function_value( &::SireFF::Point::null );
            
            Point_exposer.def( 
                "null"
                , null_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireFF::Point::operator()
        
            typedef ::SireMaths::Vector const & ( ::SireFF::Point::*__call___function_type)(  ) const;
            __call___function_type __call___function_value( &::SireFF::Point::operator() );
            
            Point_exposer.def( 
                "__call__"
                , __call___function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireFF::Point::point
        
            typedef ::SireMaths::Vector const & ( ::SireFF::Point::*point_function_type)(  ) const;
            point_function_type point_function_value( &::SireFF::Point::point );
            
            Point_exposer.def( 
                "point"
                , point_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the point in 3D space" );
        
        }
        { //::SireFF::Point::setSpace
        
            typedef void ( ::SireFF::Point::*setSpace_function_type)( ::SireVol::Space const & ) ;
            setSpace_function_type setSpace_function_value( &::SireFF::Point::setSpace );
            
            Point_exposer.def( 
                "setSpace"
                , setSpace_function_value
                , ( bp::arg("space") )
                , bp::release_gil_policy()
                , "Set the 3D space in which this point is calculated\n(although note that this 3D point, like the molecules,\nexists in the infinite cartesian space)" );
        
        }
        { //::SireFF::Point::space
        
            typedef ::SireVol::Space const & ( ::SireFF::Point::*space_function_type)(  ) const;
            space_function_type space_function_value( &::SireFF::Point::space );
            
            Point_exposer.def( 
                "space"
                , space_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the 3D space in which this point is calculated\n(although note that this 3D point, like the molecules,\nexists in the infinite cartesian space)" );
        
        }
        { //::SireFF::Point::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireFF::Point::typeName );
            
            Point_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Point::update
        
            typedef bool ( ::SireFF::Point::*update_function_type)( ::SireMol::MoleculeData const & ) ;
            update_function_type update_function_value( &::SireFF::Point::update );
            
            Point_exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("moldata") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Point::update
        
            typedef bool ( ::SireFF::Point::*update_function_type)( ::SireMol::Molecules const & ) ;
            update_function_type update_function_value( &::SireFF::Point::update );
            
            Point_exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("molecules") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Point::update
        
            typedef bool ( ::SireFF::Point::*update_function_type)( ::SireMol::MoleculeGroup const & ) ;
            update_function_type update_function_value( &::SireFF::Point::update );
            
            Point_exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("molgroup") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Point::update
        
            typedef bool ( ::SireFF::Point::*update_function_type)( ::SireMol::MolGroupsBase const & ) ;
            update_function_type update_function_value( &::SireFF::Point::update );
            
            Point_exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("molgroups") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Point::usesMoleculesIn
        
            typedef bool ( ::SireFF::Point::*usesMoleculesIn_function_type)( ::SireFF::ForceTable const & ) const;
            usesMoleculesIn_function_type usesMoleculesIn_function_value( &::SireFF::Point::usesMoleculesIn );
            
            Point_exposer.def( 
                "usesMoleculesIn"
                , usesMoleculesIn_function_value
                , ( bp::arg("forcetable") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Point::usesMoleculesIn
        
            typedef bool ( ::SireFF::Point::*usesMoleculesIn_function_type)( ::SireMol::Molecules const & ) const;
            usesMoleculesIn_function_type usesMoleculesIn_function_value( &::SireFF::Point::usesMoleculesIn );
            
            Point_exposer.def( 
                "usesMoleculesIn"
                , usesMoleculesIn_function_value
                , ( bp::arg("molecules") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Point::usesMoleculesIn
        
            typedef bool ( ::SireFF::Point::*usesMoleculesIn_function_type)( ::SireMol::MoleculeGroup const & ) const;
            usesMoleculesIn_function_type usesMoleculesIn_function_value( &::SireFF::Point::usesMoleculesIn );
            
            Point_exposer.def( 
                "usesMoleculesIn"
                , usesMoleculesIn_function_value
                , ( bp::arg("molgroup") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Point::usesMoleculesIn
        
            typedef bool ( ::SireFF::Point::*usesMoleculesIn_function_type)( ::SireMol::MolGroupsBase const & ) const;
            usesMoleculesIn_function_type usesMoleculesIn_function_value( &::SireFF::Point::usesMoleculesIn );
            
            Point_exposer.def( 
                "usesMoleculesIn"
                , usesMoleculesIn_function_value
                , ( bp::arg("molgroups") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Point::wouldUpdate
        
            typedef bool ( ::SireFF::Point::*wouldUpdate_function_type)( ::SireMol::MoleculeData const & ) const;
            wouldUpdate_function_type wouldUpdate_function_value( &::SireFF::Point::wouldUpdate );
            
            Point_exposer.def( 
                "wouldUpdate"
                , wouldUpdate_function_value
                , ( bp::arg("moldata") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Point::wouldUpdate
        
            typedef bool ( ::SireFF::Point::*wouldUpdate_function_type)( ::SireMol::Molecules const & ) const;
            wouldUpdate_function_type wouldUpdate_function_value( &::SireFF::Point::wouldUpdate );
            
            Point_exposer.def( 
                "wouldUpdate"
                , wouldUpdate_function_value
                , ( bp::arg("molecules") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Point::wouldUpdate
        
            typedef bool ( ::SireFF::Point::*wouldUpdate_function_type)( ::SireMol::MoleculeGroup const & ) const;
            wouldUpdate_function_type wouldUpdate_function_value( &::SireFF::Point::wouldUpdate );
            
            Point_exposer.def( 
                "wouldUpdate"
                , wouldUpdate_function_value
                , ( bp::arg("molgroup") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::Point::wouldUpdate
        
            typedef bool ( ::SireFF::Point::*wouldUpdate_function_type)( ::SireMol::MolGroupsBase const & ) const;
            wouldUpdate_function_type wouldUpdate_function_value( &::SireFF::Point::wouldUpdate );
            
            Point_exposer.def( 
                "wouldUpdate"
                , wouldUpdate_function_value
                , ( bp::arg("molgroups") )
                , bp::release_gil_policy()
                , "" );
        
        }
        Point_exposer.staticmethod( "areIntraMoleculePoints" );
        Point_exposer.staticmethod( "null" );
        Point_exposer.staticmethod( "typeName" );
        Point_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireFF::Point >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Point_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireFF::Point >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Point_exposer.def_pickle(sire_pickle_suite< ::SireFF::Point >());
        Point_exposer.def( "__str__", &__str__< ::SireFF::Point > );
        Point_exposer.def( "__repr__", &__str__< ::SireFF::Point > );
    }

}
