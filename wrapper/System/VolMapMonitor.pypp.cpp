// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "VolMapMonitor.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/atomcoords.h"

#include "SireMol/atomelements.h"

#include "SireMol/atommasses.h"

#include "SireMol/mgnum.h"

#include "SireMol/molecule.h"

#include "SireMol/moleculegroup.h"

#include "SireMol/partialmolecule.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireSystem/system.h"

#include "SireUnits/units.h"

#include "SireVol/gridinfo.h"

#include "volmapmonitor.h"

#include <QElapsedTimer>

#include "volmapmonitor.h"

SireSystem::VolMapMonitor __copy__(const SireSystem::VolMapMonitor &other){ return SireSystem::VolMapMonitor(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_VolMapMonitor_class(){

    { //::SireSystem::VolMapMonitor
        typedef bp::class_< SireSystem::VolMapMonitor, bp::bases< SireSystem::SystemMonitor, SireBase::Property > > VolMapMonitor_exposer_t;
        VolMapMonitor_exposer_t VolMapMonitor_exposer = VolMapMonitor_exposer_t( "VolMapMonitor", "This class create an occupation volume map showing the\nregions of space that are occupied by the monitored atoms\non average during the simulation\n\nAuthor: Christopher Woods\n", bp::init< >("Null constructor") );
        bp::scope VolMapMonitor_scope( VolMapMonitor_exposer );
        VolMapMonitor_exposer.def( bp::init< SireUnits::Dimension::Length const & >(( bp::arg("grid_spacing") ), "Construct specifying the grid spacing") );
        VolMapMonitor_exposer.def( bp::init< SireMol::MoleculeGroup const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("group"), bp::arg("map")=SireBase::PropertyMap() ), "Construct, specifying the molecule group to be monitored") );
        VolMapMonitor_exposer.def( bp::init< SireMol::MoleculeGroup const &, SireUnits::Dimension::Length const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("group"), bp::arg("grid_spacing"), bp::arg("map")=SireBase::PropertyMap() ), "Construct, specifying the molecule group to be monitored and the grid spacing\nfor the occupancy grid") );
        VolMapMonitor_exposer.def( bp::init< SireMol::MoleculeGroup const &, bool, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("group"), bp::arg("skip_light_atoms"), bp::arg("map")=SireBase::PropertyMap() ), "Construct, specifying the molecule group to be monitored\nand whether or not to ignore light atoms") );
        VolMapMonitor_exposer.def( bp::init< SireMol::MoleculeGroup const &, SireUnits::Dimension::Length const &, bool, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("group"), bp::arg("grid_spacing"), bp::arg("skip_light_atoms"), bp::arg("map")=SireBase::PropertyMap() ), "Construct, specifying the molecule group to be monitored, the grid spacing\nfor the occupancy grid, and whether or not to ignore light atoms") );
        VolMapMonitor_exposer.def( bp::init< SireSystem::VolMapMonitor const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireSystem::VolMapMonitor::averageOccupancy
        
            typedef ::QVector< float > ( ::SireSystem::VolMapMonitor::*averageOccupancy_function_type)(  ) const;
            averageOccupancy_function_type averageOccupancy_function_value( &::SireSystem::VolMapMonitor::averageOccupancy );
            
            VolMapMonitor_exposer.def( 
                "averageOccupancy"
                , averageOccupancy_function_value
                , bp::release_gil_policy()
                , "Return the average occupancy. This is a linear array that can be\naccessed using the accompanying GridInfo returned by gridInfo()" );
        
        }
        { //::SireSystem::VolMapMonitor::clearStatistics
        
            typedef void ( ::SireSystem::VolMapMonitor::*clearStatistics_function_type)(  ) ;
            clearStatistics_function_type clearStatistics_function_value( &::SireSystem::VolMapMonitor::clearStatistics );
            
            VolMapMonitor_exposer.def( 
                "clearStatistics"
                , clearStatistics_function_value
                , bp::release_gil_policy()
                , "Clear all statistics" );
        
        }
        { //::SireSystem::VolMapMonitor::gridInfo
        
            typedef ::SireVol::GridInfo ( ::SireSystem::VolMapMonitor::*gridInfo_function_type)(  ) const;
            gridInfo_function_type gridInfo_function_value( &::SireSystem::VolMapMonitor::gridInfo );
            
            VolMapMonitor_exposer.def( 
                "gridInfo"
                , gridInfo_function_value
                , bp::release_gil_policy()
                , "Return the grid dimensions" );
        
        }
        { //::SireSystem::VolMapMonitor::gridSpacing
        
            typedef ::SireUnits::Dimension::Length ( ::SireSystem::VolMapMonitor::*gridSpacing_function_type)(  ) const;
            gridSpacing_function_type gridSpacing_function_value( &::SireSystem::VolMapMonitor::gridSpacing );
            
            VolMapMonitor_exposer.def( 
                "gridSpacing"
                , gridSpacing_function_value
                , bp::release_gil_policy()
                , "Return the spacing of the grid" );
        
        }
        { //::SireSystem::VolMapMonitor::group
        
            typedef ::SireMol::MoleculeGroup const & ( ::SireSystem::VolMapMonitor::*group_function_type)(  ) const;
            group_function_type group_function_value( &::SireSystem::VolMapMonitor::group );
            
            VolMapMonitor_exposer.def( 
                "group"
                , group_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the molecule group whose atoms are being monitored" );
        
        }
        { //::SireSystem::VolMapMonitor::map
        
            typedef ::SireBase::PropertyMap ( ::SireSystem::VolMapMonitor::*map_function_type)(  ) const;
            map_function_type map_function_value( &::SireSystem::VolMapMonitor::map );
            
            VolMapMonitor_exposer.def( 
                "map"
                , map_function_value
                , bp::release_gil_policy()
                , "Return the property map used to find the properties needed by this monitor" );
        
        }
        { //::SireSystem::VolMapMonitor::monitor
        
            typedef void ( ::SireSystem::VolMapMonitor::*monitor_function_type)( ::SireSystem::System & ) ;
            monitor_function_type monitor_function_value( &::SireSystem::VolMapMonitor::monitor );
            
            VolMapMonitor_exposer.def( 
                "monitor"
                , monitor_function_value
                , ( bp::arg("system") )
                , bp::release_gil_policy()
                , "Monitor the system" );
        
        }
        { //::SireSystem::VolMapMonitor::nSamples
        
            typedef ::qint64 ( ::SireSystem::VolMapMonitor::*nSamples_function_type)(  ) const;
            nSamples_function_type nSamples_function_value( &::SireSystem::VolMapMonitor::nSamples );
            
            VolMapMonitor_exposer.def( 
                "nSamples"
                , nSamples_function_value
                , bp::release_gil_policy()
                , "Return the number of samples recorded so far" );
        
        }
        VolMapMonitor_exposer.def( bp::self != bp::self );
        { //::SireSystem::VolMapMonitor::operator=
        
            typedef ::SireSystem::VolMapMonitor & ( ::SireSystem::VolMapMonitor::*assign_function_type)( ::SireSystem::VolMapMonitor const & ) ;
            assign_function_type assign_function_value( &::SireSystem::VolMapMonitor::operator= );
            
            VolMapMonitor_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        VolMapMonitor_exposer.def( bp::self == bp::self );
        { //::SireSystem::VolMapMonitor::setGridSpacing
        
            typedef void ( ::SireSystem::VolMapMonitor::*setGridSpacing_function_type)( ::SireUnits::Dimension::Length const & ) ;
            setGridSpacing_function_type setGridSpacing_function_value( &::SireSystem::VolMapMonitor::setGridSpacing );
            
            VolMapMonitor_exposer.def( 
                "setGridSpacing"
                , setGridSpacing_function_value
                , ( bp::arg("grid_spacing") )
                , bp::release_gil_policy()
                , "Set the grid spacing on which we monitor the space" );
        
        }
        { //::SireSystem::VolMapMonitor::setGroup
        
            typedef void ( ::SireSystem::VolMapMonitor::*setGroup_function_type)( ::SireMol::MoleculeGroup const &,::SireBase::PropertyMap const & ) ;
            setGroup_function_type setGroup_function_value( &::SireSystem::VolMapMonitor::setGroup );
            
            VolMapMonitor_exposer.def( 
                "setGroup"
                , setGroup_function_value
                , ( bp::arg("group"), bp::arg("map")=SireBase::PropertyMap() )
                , "Set the molecule group being monitored (together with the property map\nused to get the required properties from that group)" );
        
        }
        { //::SireSystem::VolMapMonitor::setPropertyMap
        
            typedef void ( ::SireSystem::VolMapMonitor::*setPropertyMap_function_type)( ::SireBase::PropertyMap const & ) ;
            setPropertyMap_function_type setPropertyMap_function_value( &::SireSystem::VolMapMonitor::setPropertyMap );
            
            VolMapMonitor_exposer.def( 
                "setPropertyMap"
                , setPropertyMap_function_value
                , ( bp::arg("map") )
                , bp::release_gil_policy()
                , "Set the property map to be used by this monitor" );
        
        }
        { //::SireSystem::VolMapMonitor::setSkipLightAtoms
        
            typedef void ( ::SireSystem::VolMapMonitor::*setSkipLightAtoms_function_type)( bool ) ;
            setSkipLightAtoms_function_type setSkipLightAtoms_function_value( &::SireSystem::VolMapMonitor::setSkipLightAtoms );
            
            VolMapMonitor_exposer.def( 
                "setSkipLightAtoms"
                , setSkipLightAtoms_function_value
                , ( bp::arg("on") )
                , bp::release_gil_policy()
                , "Turn on skipping light atoms" );
        
        }
        { //::SireSystem::VolMapMonitor::skippingLightAtoms
        
            typedef bool ( ::SireSystem::VolMapMonitor::*skippingLightAtoms_function_type)(  ) const;
            skippingLightAtoms_function_type skippingLightAtoms_function_value( &::SireSystem::VolMapMonitor::skippingLightAtoms );
            
            VolMapMonitor_exposer.def( 
                "skippingLightAtoms"
                , skippingLightAtoms_function_value
                , bp::release_gil_policy()
                , "Return whether or not we are skipping light atoms (e.g. hydrogen)" );
        
        }
        { //::SireSystem::VolMapMonitor::toString
        
            typedef ::QString ( ::SireSystem::VolMapMonitor::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireSystem::VolMapMonitor::toString );
            
            VolMapMonitor_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireSystem::VolMapMonitor::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireSystem::VolMapMonitor::typeName );
            
            VolMapMonitor_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireSystem::VolMapMonitor::volumeMap
        
            typedef ::SireMol::VolumeMap ( ::SireSystem::VolMapMonitor::*volumeMap_function_type)(  ) const;
            volumeMap_function_type volumeMap_function_value( &::SireSystem::VolMapMonitor::volumeMap );
            
            VolMapMonitor_exposer.def( 
                "volumeMap"
                , volumeMap_function_value
                , bp::release_gil_policy()
                , "Return the actual volume map" );
        
        }
        { //::SireSystem::VolMapMonitor::what
        
            typedef char const * ( ::SireSystem::VolMapMonitor::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireSystem::VolMapMonitor::what );
            
            VolMapMonitor_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        VolMapMonitor_exposer.staticmethod( "typeName" );
        VolMapMonitor_exposer.def( "__copy__", &__copy__);
        VolMapMonitor_exposer.def( "__deepcopy__", &__copy__);
        VolMapMonitor_exposer.def( "clone", &__copy__);
        VolMapMonitor_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireSystem::VolMapMonitor >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        VolMapMonitor_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireSystem::VolMapMonitor >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        VolMapMonitor_exposer.def_pickle(sire_pickle_suite< ::SireSystem::VolMapMonitor >());
        VolMapMonitor_exposer.def( "__str__", &__str__< ::SireSystem::VolMapMonitor > );
        VolMapMonitor_exposer.def( "__repr__", &__str__< ::SireSystem::VolMapMonitor > );
    }

}
