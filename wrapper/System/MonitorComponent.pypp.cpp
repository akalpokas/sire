// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "MonitorComponent.pypp.hpp"

namespace bp = boost::python;

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "monitorcomponent.h"

#include "system.h"

#include "monitorcomponent.h"

SireSystem::MonitorComponent __copy__(const SireSystem::MonitorComponent &other){ return SireSystem::MonitorComponent(other); }

#include "Helpers/copy.hpp"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_MonitorComponent_class(){

    { //::SireSystem::MonitorComponent
        typedef bp::class_< SireSystem::MonitorComponent, bp::bases< SireSystem::SystemMonitor, SireBase::Property > > MonitorComponent_exposer_t;
        MonitorComponent_exposer_t MonitorComponent_exposer = MonitorComponent_exposer_t( "MonitorComponent", "This monitor is used to monitor the value of a component of the system.\nIt can be used to generate an average of that component. It could be\nused, for example, to generate average energies, or accumulate\nfree energies during a simulation\n\nAuthor: Christopher Woods\n", bp::init< >("Null constructor") );
        bp::scope MonitorComponent_scope( MonitorComponent_exposer );
        MonitorComponent_exposer.def( bp::init< SireCAS::Symbol const & >(( bp::arg("component") ), "Construct a monitor to collect the average value of the component\nrepresented by the symbol component") );
        MonitorComponent_exposer.def( bp::init< SireCAS::Symbol const &, SireMaths::Accumulator const & >(( bp::arg("component"), bp::arg("accumulator") ), "Construct a monitor to accumulate the value of the component\nrepresented by the symbol component using the accumulator\nin accumulator") );
        MonitorComponent_exposer.def( bp::init< SireSystem::MonitorComponent const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireSystem::MonitorComponent::accumulator
        
            typedef ::SireMaths::Accumulator const & ( ::SireSystem::MonitorComponent::*accumulator_function_type)(  ) const;
            accumulator_function_type accumulator_function_value( &::SireSystem::MonitorComponent::accumulator );
            
            MonitorComponent_exposer.def( 
                "accumulator"
                , accumulator_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the accumulator that is being used to accumulate the\nvalues of the component being monitored" );
        
        }
        { //::SireSystem::MonitorComponent::clearStatistics
        
            typedef void ( ::SireSystem::MonitorComponent::*clearStatistics_function_type)(  ) ;
            clearStatistics_function_type clearStatistics_function_value( &::SireSystem::MonitorComponent::clearStatistics );
            
            MonitorComponent_exposer.def( 
                "clearStatistics"
                , clearStatistics_function_value
                , bp::release_gil_policy()
                , "Clear the statistics in this monitor" );
        
        }
        { //::SireSystem::MonitorComponent::component
        
            typedef ::SireCAS::Symbol const & ( ::SireSystem::MonitorComponent::*component_function_type)(  ) const;
            component_function_type component_function_value( &::SireSystem::MonitorComponent::component );
            
            MonitorComponent_exposer.def( 
                "component"
                , component_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the symbol representing the component being monitored" );
        
        }
        { //::SireSystem::MonitorComponent::monitor
        
            typedef void ( ::SireSystem::MonitorComponent::*monitor_function_type)( ::SireSystem::System & ) ;
            monitor_function_type monitor_function_value( &::SireSystem::MonitorComponent::monitor );
            
            MonitorComponent_exposer.def( 
                "monitor"
                , monitor_function_value
                , ( bp::arg("system") )
                , bp::release_gil_policy()
                , "Call this function to add the statistics of the monitored\ncomponent from the passed system to the accumulator" );
        
        }
        MonitorComponent_exposer.def( bp::self != bp::self );
        { //::SireSystem::MonitorComponent::operator=
        
            typedef ::SireSystem::MonitorComponent & ( ::SireSystem::MonitorComponent::*assign_function_type)( ::SireSystem::MonitorComponent const & ) ;
            assign_function_type assign_function_value( &::SireSystem::MonitorComponent::operator= );
            
            MonitorComponent_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        MonitorComponent_exposer.def( bp::self == bp::self );
        { //::SireSystem::MonitorComponent::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireSystem::MonitorComponent::typeName );
            
            MonitorComponent_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        MonitorComponent_exposer.staticmethod( "typeName" );
        MonitorComponent_exposer.def( "__copy__", &__copy__<SireSystem::MonitorComponent>);
        MonitorComponent_exposer.def( "__deepcopy__", &__copy__<SireSystem::MonitorComponent>);
        MonitorComponent_exposer.def( "clone", &__copy__<SireSystem::MonitorComponent>);
        MonitorComponent_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireSystem::MonitorComponent >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        MonitorComponent_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireSystem::MonitorComponent >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        MonitorComponent_exposer.def_pickle(sire_pickle_suite< ::SireSystem::MonitorComponent >());
        MonitorComponent_exposer.def( "__str__", &__str__< ::SireSystem::MonitorComponent > );
        MonitorComponent_exposer.def( "__repr__", &__str__< ::SireSystem::MonitorComponent > );
    }

}
