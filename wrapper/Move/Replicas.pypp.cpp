// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "Replicas.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireID/index.h"

#include "SireMaths/rangenerator.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "replica.h"

#include "replicas.h"

#include "replicas.h"

SireMove::Replicas __copy__(const SireMove::Replicas &other){ return SireMove::Replicas(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_Replicas_class(){

    { //::SireMove::Replicas
        typedef bp::class_< SireMove::Replicas, bp::bases< SireMove::SupraSystem, SireBase::Property > > Replicas_exposer_t;
        Replicas_exposer_t Replicas_exposer = Replicas_exposer_t( "Replicas", "This class is used to hold all of the replicas in\na replica exchange simulation\n\nAuthor: Christopher Woods\n", bp::init< >("Constructor") );
        bp::scope Replicas_scope( Replicas_exposer );
        Replicas_exposer.def( bp::init< int >(( bp::arg("n") ), "Construct a set of n replicas") );
        Replicas_exposer.def( bp::init< SireSystem::System const &, bp::optional< int > >(( bp::arg("system"), bp::arg("n")=(int)(1) ), "Construct a set of n replicas that contain the sub-system system") );
        Replicas_exposer.def( bp::init< QVector< SireSystem::System > const & >(( bp::arg("systems") ), "Construct a set of replicas from the passed array of systems") );
        Replicas_exposer.def( bp::init< SireMove::SupraSubSystem const &, bp::optional< int > >(( bp::arg("subsystem"), bp::arg("n")=(int)(1) ), "Construct a Replicas object that contains n copies of the\npassed subsystem") );
        Replicas_exposer.def( bp::init< SireMove::SupraSystem const & >(( bp::arg("suprasystem") ), "Construct from the passed SupraSystem") );
        Replicas_exposer.def( bp::init< SireMove::Replicas const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMove::Replicas::at
        
            typedef ::SireMove::Replica const & ( ::SireMove::Replicas::*at_function_type)( int ) const;
            at_function_type at_function_value( &::SireMove::Replicas::at );
            
            Replicas_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the ith replica\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMove::Replicas::collectSupraStats
        
            typedef void ( ::SireMove::Replicas::*collectSupraStats_function_type)(  ) ;
            collectSupraStats_function_type collectSupraStats_function_value( &::SireMove::Replicas::collectSupraStats );
            
            Replicas_exposer.def( 
                "collectSupraStats"
                , collectSupraStats_function_value
                , bp::release_gil_policy()
                , "Collect statistics - this records the current lambdaTrajectory and adds it to the history" );
        
        }
        { //::SireMove::Replicas::lambdaTrajectory
        
            typedef ::QVector< double > ( ::SireMove::Replicas::*lambdaTrajectory_function_type)(  ) const;
            lambdaTrajectory_function_type lambdaTrajectory_function_value( &::SireMove::Replicas::lambdaTrajectory );
            
            Replicas_exposer.def( 
                "lambdaTrajectory"
                , lambdaTrajectory_function_value
                , bp::release_gil_policy()
                , "Return the lambda values for each of the replicas, in replica ID order.\nThis allows the lambda trajectory for each replica to be easily\ncollected during a simulation" );
        
        }
        { //::SireMove::Replicas::lambdaTrajectoryHistory
        
            typedef ::QList< QVector< double > > ( ::SireMove::Replicas::*lambdaTrajectoryHistory_function_type)(  ) const;
            lambdaTrajectoryHistory_function_type lambdaTrajectoryHistory_function_value( &::SireMove::Replicas::lambdaTrajectoryHistory );
            
            Replicas_exposer.def( 
                "lambdaTrajectoryHistory"
                , lambdaTrajectoryHistory_function_value
                , bp::release_gil_policy()
                , "Return the history of lambda values sampled by each replica" );
        
        }
        { //::SireMove::Replicas::nReplicas
        
            typedef int ( ::SireMove::Replicas::*nReplicas_function_type)(  ) const;
            nReplicas_function_type nReplicas_function_value( &::SireMove::Replicas::nReplicas );
            
            Replicas_exposer.def( 
                "nReplicas"
                , nReplicas_function_value
                , bp::release_gil_policy()
                , "Return the number of replicas in this set" );
        
        }
        Replicas_exposer.def( bp::self != bp::self );
        { //::SireMove::Replicas::operator=
        
            typedef ::SireMove::Replicas & ( ::SireMove::Replicas::*assign_function_type)( ::SireMove::Replicas const & ) ;
            assign_function_type assign_function_value( &::SireMove::Replicas::operator= );
            
            Replicas_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        Replicas_exposer.def( bp::self == bp::self );
        { //::SireMove::Replicas::operator[]
        
            typedef ::SireMove::Replica const & ( ::SireMove::Replicas::*__getitem___function_type)( int ) const;
            __getitem___function_type __getitem___function_value( &::SireMove::Replicas::operator[] );
            
            Replicas_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMove::Replicas::replicaIDs
        
            typedef ::QVector< unsigned int > const & ( ::SireMove::Replicas::*replicaIDs_function_type)(  ) const;
            replicaIDs_function_type replicaIDs_function_value( &::SireMove::Replicas::replicaIDs );
            
            Replicas_exposer.def( 
                "replicaIDs"
                , replicaIDs_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the current IDs of the replicas (in order of replica index)" );
        
        }
        { //::SireMove::Replicas::resetReplicaIDs
        
            typedef void ( ::SireMove::Replicas::*resetReplicaIDs_function_type)(  ) ;
            resetReplicaIDs_function_type resetReplicaIDs_function_value( &::SireMove::Replicas::resetReplicaIDs );
            
            Replicas_exposer.def( 
                "resetReplicaIDs"
                , resetReplicaIDs_function_value
                , bp::release_gil_policy()
                , "Reset the replica IDs - this sets the ID of the ith replica to i" );
        
        }
        { //::SireMove::Replicas::setChemicalPotential
        
            typedef void ( ::SireMove::Replicas::*setChemicalPotential_function_type)( ::SireUnits::Dimension::MolarEnergy const & ) ;
            setChemicalPotential_function_type setChemicalPotential_function_value( &::SireMove::Replicas::setChemicalPotential );
            
            Replicas_exposer.def( 
                "setChemicalPotential"
                , setChemicalPotential_function_value
                , ( bp::arg("chemical_potential") )
                , bp::release_gil_policy()
                , "Set the chemical potential of the ensemble sampled by all\nreplicas to chemical_potential\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireMove::Replicas::setChemicalPotential
        
            typedef void ( ::SireMove::Replicas::*setChemicalPotential_function_type)( int,::SireUnits::Dimension::MolarEnergy const & ) ;
            setChemicalPotential_function_type setChemicalPotential_function_value( &::SireMove::Replicas::setChemicalPotential );
            
            Replicas_exposer.def( 
                "setChemicalPotential"
                , setChemicalPotential_function_value
                , ( bp::arg("i"), bp::arg("chemical_potential") )
                , bp::release_gil_policy()
                , "Set the chemical potential of the ensemble sampled by the ith\nreplica to chemical_potential\nThrow: SireError::invalid_index\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireMove::Replicas::setEnergyComponent
        
            typedef void ( ::SireMove::Replicas::*setEnergyComponent_function_type)( ::SireCAS::Symbol const & ) ;
            setEnergyComponent_function_type setEnergyComponent_function_value( &::SireMove::Replicas::setEnergyComponent );
            
            Replicas_exposer.def( 
                "setEnergyComponent"
                , setEnergyComponent_function_value
                , ( bp::arg("symbol") )
                , bp::release_gil_policy()
                , "Set the energy component sampled by all replicas to symbol" );
        
        }
        { //::SireMove::Replicas::setEnergyComponent
        
            typedef void ( ::SireMove::Replicas::*setEnergyComponent_function_type)( int,::SireCAS::Symbol const & ) ;
            setEnergyComponent_function_type setEnergyComponent_function_value( &::SireMove::Replicas::setEnergyComponent );
            
            Replicas_exposer.def( 
                "setEnergyComponent"
                , setEnergyComponent_function_value
                , ( bp::arg("i"), bp::arg("symbol") )
                , bp::release_gil_policy()
                , "Set the energy component sampled by the ith replica to symbol\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMove::Replicas::setFugacity
        
            typedef void ( ::SireMove::Replicas::*setFugacity_function_type)( ::SireUnits::Dimension::Pressure const & ) ;
            setFugacity_function_type setFugacity_function_value( &::SireMove::Replicas::setFugacity );
            
            Replicas_exposer.def( 
                "setFugacity"
                , setFugacity_function_value
                , ( bp::arg("fugacity") )
                , bp::release_gil_policy()
                , "Set the fugacity of the ensemble sampled by all replicas\nto fugacity\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireMove::Replicas::setFugacity
        
            typedef void ( ::SireMove::Replicas::*setFugacity_function_type)( int,::SireUnits::Dimension::Pressure const & ) ;
            setFugacity_function_type setFugacity_function_value( &::SireMove::Replicas::setFugacity );
            
            Replicas_exposer.def( 
                "setFugacity"
                , setFugacity_function_value
                , ( bp::arg("i"), bp::arg("fugacity") )
                , bp::release_gil_policy()
                , "Set the fugacity of the ensemble sampled by the ith replica\nto fugacity\nThrow: SireError::invalid_index\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireMove::Replicas::setGenerator
        
            typedef void ( ::SireMove::Replicas::*setGenerator_function_type)( ::SireMaths::RanGenerator const & ) ;
            setGenerator_function_type setGenerator_function_value( &::SireMove::Replicas::setGenerator );
            
            Replicas_exposer.def( 
                "setGenerator"
                , setGenerator_function_value
                , ( bp::arg("rangenerator") )
                , bp::release_gil_policy()
                , "Set the random number generator used by all of the replicas\nto rangenerator - this doesnt give all of the replicas\nthe same generator - rather it uses this generator to\nreproducibly generate new generators for each replica" );
        
        }
        { //::SireMove::Replicas::setGenerator
        
            typedef void ( ::SireMove::Replicas::*setGenerator_function_type)( int,::SireMaths::RanGenerator const & ) ;
            setGenerator_function_type setGenerator_function_value( &::SireMove::Replicas::setGenerator );
            
            Replicas_exposer.def( 
                "setGenerator"
                , setGenerator_function_value
                , ( bp::arg("i"), bp::arg("rangenerator") )
                , bp::release_gil_policy()
                , "Set the random number generator used by the ith replica to generator\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMove::Replicas::setLambdaComponent
        
            typedef void ( ::SireMove::Replicas::*setLambdaComponent_function_type)( ::SireCAS::Symbol const & ) ;
            setLambdaComponent_function_type setLambdaComponent_function_value( &::SireMove::Replicas::setLambdaComponent );
            
            Replicas_exposer.def( 
                "setLambdaComponent"
                , setLambdaComponent_function_value
                , ( bp::arg("symbol") )
                , bp::release_gil_policy()
                , "Set the lambda component used for Hamiltonian replica exchange\nfor all replicas to symbol" );
        
        }
        { //::SireMove::Replicas::setLambdaComponent
        
            typedef void ( ::SireMove::Replicas::*setLambdaComponent_function_type)( int,::SireCAS::Symbol const & ) ;
            setLambdaComponent_function_type setLambdaComponent_function_value( &::SireMove::Replicas::setLambdaComponent );
            
            Replicas_exposer.def( 
                "setLambdaComponent"
                , setLambdaComponent_function_value
                , ( bp::arg("i"), bp::arg("symbol") )
                , bp::release_gil_policy()
                , "Set the lambda component used by the ith replica for Hamiltonian\nreplica exchange to symbol\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMove::Replicas::setLambdaValue
        
            typedef void ( ::SireMove::Replicas::*setLambdaValue_function_type)( double ) ;
            setLambdaValue_function_type setLambdaValue_function_value( &::SireMove::Replicas::setLambdaValue );
            
            Replicas_exposer.def( 
                "setLambdaValue"
                , setLambdaValue_function_value
                , ( bp::arg("value") )
                , bp::release_gil_policy()
                , "Set the lambda value for Hamiltonian replica exchange for all\nreplicas to value" );
        
        }
        { //::SireMove::Replicas::setLambdaValue
        
            typedef void ( ::SireMove::Replicas::*setLambdaValue_function_type)( int,double ) ;
            setLambdaValue_function_type setLambdaValue_function_value( &::SireMove::Replicas::setLambdaValue );
            
            Replicas_exposer.def( 
                "setLambdaValue"
                , setLambdaValue_function_value
                , ( bp::arg("i"), bp::arg("value") )
                , bp::release_gil_policy()
                , "Set the value of lambda used for Hamiltonian replica exchange for\nthe ith replica to value\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMove::Replicas::setPressure
        
            typedef void ( ::SireMove::Replicas::*setPressure_function_type)( ::SireUnits::Dimension::Pressure const & ) ;
            setPressure_function_type setPressure_function_value( &::SireMove::Replicas::setPressure );
            
            Replicas_exposer.def( 
                "setPressure"
                , setPressure_function_value
                , ( bp::arg("pressure") )
                , bp::release_gil_policy()
                , "Set the pressure of the ensemble sampled by all replicas\nto pressure" );
        
        }
        { //::SireMove::Replicas::setPressure
        
            typedef void ( ::SireMove::Replicas::*setPressure_function_type)( int,::SireUnits::Dimension::Pressure const & ) ;
            setPressure_function_type setPressure_function_value( &::SireMove::Replicas::setPressure );
            
            Replicas_exposer.def( 
                "setPressure"
                , setPressure_function_value
                , ( bp::arg("i"), bp::arg("pressure") )
                , bp::release_gil_policy()
                , "Set the pressure of the ensemble sampled by the ith replica\nto pressure\nThrow: SireError::invalid_index\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireMove::Replicas::setReplica
        
            typedef void ( ::SireMove::Replicas::*setReplica_function_type)( ::SireMove::Replica const & ) ;
            setReplica_function_type setReplica_function_value( &::SireMove::Replicas::setReplica );
            
            Replicas_exposer.def( 
                "setReplica"
                , setReplica_function_value
                , ( bp::arg("replica") )
                , bp::release_gil_policy()
                , "Set the replicas to all be a copy of replica" );
        
        }
        { //::SireMove::Replicas::setReplica
        
            typedef void ( ::SireMove::Replicas::*setReplica_function_type)( int,::SireMove::Replica const & ) ;
            setReplica_function_type setReplica_function_value( &::SireMove::Replicas::setReplica );
            
            Replicas_exposer.def( 
                "setReplica"
                , setReplica_function_value
                , ( bp::arg("i"), bp::arg("replica") )
                , bp::release_gil_policy()
                , "Set the ith replica equal to replica\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMove::Replicas::setReplicas
        
            typedef void ( ::SireMove::Replicas::*setReplicas_function_type)( ::SireMove::Replicas const & ) ;
            setReplicas_function_type setReplicas_function_value( &::SireMove::Replicas::setReplicas );
            
            Replicas_exposer.def( 
                "setReplicas"
                , setReplicas_function_value
                , ( bp::arg("replicas") )
                , bp::release_gil_policy()
                , "Set the replicas from a copy of passed replicas" );
        
        }
        { //::SireMove::Replicas::setSpaceProperty
        
            typedef void ( ::SireMove::Replicas::*setSpaceProperty_function_type)( ::SireBase::PropertyName const & ) ;
            setSpaceProperty_function_type setSpaceProperty_function_value( &::SireMove::Replicas::setSpaceProperty );
            
            Replicas_exposer.def( 
                "setSpaceProperty"
                , setSpaceProperty_function_value
                , ( bp::arg("spaceproperty") )
                , bp::release_gil_policy()
                , "Set the property used to find the simulation space for all replicas\nto spaceproperty" );
        
        }
        { //::SireMove::Replicas::setSpaceProperty
        
            typedef void ( ::SireMove::Replicas::*setSpaceProperty_function_type)( int,::SireBase::PropertyName const & ) ;
            setSpaceProperty_function_type setSpaceProperty_function_value( &::SireMove::Replicas::setSpaceProperty );
            
            Replicas_exposer.def( 
                "setSpaceProperty"
                , setSpaceProperty_function_value
                , ( bp::arg("i"), bp::arg("spaceproperty") )
                , bp::release_gil_policy()
                , "Set the space property used by the ith replica to spaceproperty\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMove::Replicas::setSubSystem
        
            typedef void ( ::SireMove::Replicas::*setSubSystem_function_type)( ::SireSystem::System const & ) ;
            setSubSystem_function_type setSubSystem_function_value( &::SireMove::Replicas::setSubSystem );
            
            Replicas_exposer.def( 
                "setSubSystem"
                , setSubSystem_function_value
                , ( bp::arg("system") )
                , bp::release_gil_policy()
                , "Overloaded function used to set all of the replicas equal to subsystem" );
        
        }
        { //::SireMove::Replicas::setSubSystem
        
            typedef void ( ::SireMove::Replicas::*setSubSystem_function_type)( ::SireMove::SupraSubSystem const & ) ;
            setSubSystem_function_type setSubSystem_function_value( &::SireMove::Replicas::setSubSystem );
            
            Replicas_exposer.def( 
                "setSubSystem"
                , setSubSystem_function_value
                , ( bp::arg("subsystem") )
                , bp::release_gil_policy()
                , "Overloaded function used to set all of the replicas equal to subsystem" );
        
        }
        { //::SireMove::Replicas::setSubSystem
        
            typedef void ( ::SireMove::Replicas::*setSubSystem_function_type)( int,::SireSystem::System const & ) ;
            setSubSystem_function_type setSubSystem_function_value( &::SireMove::Replicas::setSubSystem );
            
            Replicas_exposer.def( 
                "setSubSystem"
                , setSubSystem_function_value
                , ( bp::arg("i"), bp::arg("system") )
                , bp::release_gil_policy()
                , "Overloaded function used to set the ith subsystem equal to system\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMove::Replicas::setSubSystem
        
            typedef void ( ::SireMove::Replicas::*setSubSystem_function_type)( int,::SireMove::SupraSubSystem const & ) ;
            setSubSystem_function_type setSubSystem_function_value( &::SireMove::Replicas::setSubSystem );
            
            Replicas_exposer.def( 
                "setSubSystem"
                , setSubSystem_function_value
                , ( bp::arg("i"), bp::arg("subsystem") )
                , bp::release_gil_policy()
                , "Overloaded function used to ensure that this system only contains\nReplica subsystems\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMove::Replicas::setTemperature
        
            typedef void ( ::SireMove::Replicas::*setTemperature_function_type)( ::SireUnits::Dimension::Temperature const & ) ;
            setTemperature_function_type setTemperature_function_value( &::SireMove::Replicas::setTemperature );
            
            Replicas_exposer.def( 
                "setTemperature"
                , setTemperature_function_value
                , ( bp::arg("temperature") )
                , bp::release_gil_policy()
                , "Set the temperature of the ensemble sampled by all replicas\nto temperature\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireMove::Replicas::setTemperature
        
            typedef void ( ::SireMove::Replicas::*setTemperature_function_type)( int,::SireUnits::Dimension::Temperature const & ) ;
            setTemperature_function_type setTemperature_function_value( &::SireMove::Replicas::setTemperature );
            
            Replicas_exposer.def( 
                "setTemperature"
                , setTemperature_function_value
                , ( bp::arg("i"), bp::arg("temperature") )
                , bp::release_gil_policy()
                , "Set the temperature of the ensemble sampled by\nthe ith replica to temperature\nThrow: SireError::invalid_index\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireMove::Replicas::swapMolecules
        
            typedef void ( ::SireMove::Replicas::*swapMolecules_function_type)( int,int ) ;
            swapMolecules_function_type swapMolecules_function_value( &::SireMove::Replicas::swapMolecules );
            
            Replicas_exposer.def( 
                "swapMolecules"
                , swapMolecules_function_value
                , ( bp::arg("i"), bp::arg("j") )
                , bp::release_gil_policy()
                , "Swap the molecules between replicas i and j\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMove::Replicas::swapSystems
        
            typedef void ( ::SireMove::Replicas::*swapSystems_function_type)( int,int,bool ) ;
            swapSystems_function_type swapSystems_function_value( &::SireMove::Replicas::swapSystems );
            
            Replicas_exposer.def( 
                "swapSystems"
                , swapSystems_function_value
                , ( bp::arg("i"), bp::arg("j"), bp::arg("swap_monitors")=(bool)(true) )
                , "Swap the systems between replicas i and j. If swap_monitors is\ntrue then the monitors are swapped as well\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMove::Replicas::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMove::Replicas::typeName );
            
            Replicas_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        Replicas_exposer.staticmethod( "typeName" );
        Replicas_exposer.def( "__copy__", &__copy__);
        Replicas_exposer.def( "__deepcopy__", &__copy__);
        Replicas_exposer.def( "clone", &__copy__);
        Replicas_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMove::Replicas >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Replicas_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMove::Replicas >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Replicas_exposer.def_pickle(sire_pickle_suite< ::SireMove::Replicas >());
        Replicas_exposer.def( "__str__", &__str__< ::SireMove::Replicas > );
        Replicas_exposer.def( "__repr__", &__str__< ::SireMove::Replicas > );
        Replicas_exposer.def( "__len__", &__len_size< ::SireMove::Replicas > );
    }

}
