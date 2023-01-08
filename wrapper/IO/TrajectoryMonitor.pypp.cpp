// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "TrajectoryMonitor.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireMol/moleculegroup.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireSystem/system.h"

#include "pdb.h"

#include "trajectorymonitor.h"

#include <QByteArray>

#include <QDebug>

#include <QDir>

#include <QFile>

#include <QFileInfo>

#include <QTemporaryFile>

#include "trajectorymonitor.h"

SireIO::TrajectoryMonitor __copy__(const SireIO::TrajectoryMonitor &other){ return SireIO::TrajectoryMonitor(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_TrajectoryMonitor_class(){

    { //::SireIO::TrajectoryMonitor
        typedef bp::class_< SireIO::TrajectoryMonitor, bp::bases< SireSystem::SystemMonitor, SireBase::Property > > TrajectoryMonitor_exposer_t;
        TrajectoryMonitor_exposer_t TrajectoryMonitor_exposer = TrajectoryMonitor_exposer_t( "TrajectoryMonitor", "This is a monitor that can be used to save a trajectory\nof an arbitrary collection of molecules from the system\n\nAuthor: Christopher Woods\n", bp::init< >("Null constructor") );
        bp::scope TrajectoryMonitor_scope( TrajectoryMonitor_exposer );
        TrajectoryMonitor_exposer.def( bp::init< SireMol::MoleculeGroup const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("molgroup"), bp::arg("map")=SireBase::PropertyMap() ), "Construct a monitor that monitors the trajectory of the molecules\nin the molecule group molgroup. This writes the trajectory using\nthe PDB writer, and uses the (optionally) supplied property map to\ncontrol what is written") );
        TrajectoryMonitor_exposer.def( bp::init< SireMol::MoleculeGroup const &, SireIO::IOBase const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("molgroup"), bp::arg("writer"), bp::arg("map")=SireBase::PropertyMap() ), "Construct a monitor that monitors the trajectory of the molecules in\nthe molecule group molgroup, writing the trajectory using the\nmolecule write writer, and using the (optionally) supplied property\nmap to control what is written") );
        TrajectoryMonitor_exposer.def( bp::init< SireMol::MGID const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("mgid"), bp::arg("map")=SireBase::PropertyMap() ), "Construct a monitor that monitors the trajectory of the molecules\nin the molecule group with ID mgid. This writes the trajectory using\nthe PDB writer, and uses the (optionally) supplied property map to\ncontrol what is written") );
        TrajectoryMonitor_exposer.def( bp::init< SireMol::MGID const &, SireIO::IOBase const &, bp::optional< SireBase::PropertyMap const & > >(( bp::arg("mgid"), bp::arg("writer"), bp::arg("map")=SireBase::PropertyMap() ), "Construct a monitor that monitors the trajectory of the molecules in\nthe molecule group with ID mgid, writing the trajectory using the\nmolecule write writer, and using the (optionally) supplied property\nmap to control what is written") );
        TrajectoryMonitor_exposer.def( bp::init< SireIO::TrajectoryMonitor const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireIO::TrajectoryMonitor::clearStatistics
        
            typedef void ( ::SireIO::TrajectoryMonitor::*clearStatistics_function_type)(  ) ;
            clearStatistics_function_type clearStatistics_function_value( &::SireIO::TrajectoryMonitor::clearStatistics );
            
            TrajectoryMonitor_exposer.def( 
                "clearStatistics"
                , clearStatistics_function_value
                , bp::release_gil_policy()
                , "Clear all statistics - this will clear all frames of the trajectory" );
        
        }
        { //::SireIO::TrajectoryMonitor::monitor
        
            typedef void ( ::SireIO::TrajectoryMonitor::*monitor_function_type)( ::SireSystem::System & ) ;
            monitor_function_type monitor_function_value( &::SireIO::TrajectoryMonitor::monitor );
            
            TrajectoryMonitor_exposer.def( 
                "monitor"
                , monitor_function_value
                , ( bp::arg("system") )
                , bp::release_gil_policy()
                , "Monitor the system, writing an additional frame of the trajectory\nto this monitor" );
        
        }
        TrajectoryMonitor_exposer.def( bp::self != bp::self );
        { //::SireIO::TrajectoryMonitor::operator=
        
            typedef ::SireIO::TrajectoryMonitor & ( ::SireIO::TrajectoryMonitor::*assign_function_type)( ::SireIO::TrajectoryMonitor const & ) ;
            assign_function_type assign_function_value( &::SireIO::TrajectoryMonitor::operator= );
            
            TrajectoryMonitor_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        TrajectoryMonitor_exposer.def( bp::self == bp::self );
        { //::SireIO::TrajectoryMonitor::setTempDir
        
            typedef void ( ::SireIO::TrajectoryMonitor::*setTempDir_function_type)( ::QString const & ) ;
            setTempDir_function_type setTempDir_function_value( &::SireIO::TrajectoryMonitor::setTempDir );
            
            TrajectoryMonitor_exposer.def( 
                "setTempDir"
                , setTempDir_function_value
                , ( bp::arg("tempdir") )
                , bp::release_gil_policy()
                , "Set the temporary directory used to store the trajectory as it\nis being monitored during the simulation" );
        
        }
        { //::SireIO::TrajectoryMonitor::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireIO::TrajectoryMonitor::typeName );
            
            TrajectoryMonitor_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireIO::TrajectoryMonitor::writeToDisk
        
            typedef void ( ::SireIO::TrajectoryMonitor::*writeToDisk_function_type)( ::QString const & ) const;
            writeToDisk_function_type writeToDisk_function_value( &::SireIO::TrajectoryMonitor::writeToDisk );
            
            TrajectoryMonitor_exposer.def( 
                "writeToDisk"
                , writeToDisk_function_value
                , ( bp::arg("file_template") )
                , bp::release_gil_policy()
                , "Write the trajectory to disk, using the file template file_template.\nThis writes each frame to a separate file, numbering them sequentially\nfrom 0, replacing XXXXXX with the frame number" );
        
        }
        TrajectoryMonitor_exposer.staticmethod( "typeName" );
        TrajectoryMonitor_exposer.def( "__copy__", &__copy__);
        TrajectoryMonitor_exposer.def( "__deepcopy__", &__copy__);
        TrajectoryMonitor_exposer.def( "clone", &__copy__);
        TrajectoryMonitor_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireIO::TrajectoryMonitor >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        TrajectoryMonitor_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireIO::TrajectoryMonitor >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        TrajectoryMonitor_exposer.def_pickle(sire_pickle_suite< ::SireIO::TrajectoryMonitor >());
        TrajectoryMonitor_exposer.def( "__str__", &__str__< ::SireIO::TrajectoryMonitor > );
        TrajectoryMonitor_exposer.def( "__repr__", &__str__< ::SireIO::TrajectoryMonitor > );
    }

}
