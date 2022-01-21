// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "Replica.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireFF/errors.h"

#include "SireMaths/rangenerator.h"

#include "SireMol/molecules.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireSystem/system.h"

#include "SireUnits/units.h"

#include "replica.h"

#include "replicas.h"

#include "replica.h"

SireMove::Replica __copy__(const SireMove::Replica &other){ return SireMove::Replica(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_Replica_class(){

    { //::SireMove::Replica
        typedef bp::class_< SireMove::Replica, bp::bases< SireMove::SupraSubSystem, SireBase::Property > > Replica_exposer_t;
        Replica_exposer_t Replica_exposer = Replica_exposer_t( "Replica", "This is a replica within a replica exchange simulation\n\nAuthor: Christopher Woods\n", bp::init< >("Constructor") );
        bp::scope Replica_scope( Replica_exposer );
        Replica_exposer.def( bp::init< SireMove::SupraSubSystem const & >(( bp::arg("subsystem") ), "Construct from another SubSystem") );
        Replica_exposer.def( bp::init< SireMove::Replica const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMove::Replica::chemicalPotential
        
            typedef ::SireUnits::Dimension::MolarEnergy ( ::SireMove::Replica::*chemicalPotential_function_type)(  ) const;
            chemicalPotential_function_type chemicalPotential_function_value( &::SireMove::Replica::chemicalPotential );
            
            Replica_exposer.def( 
                "chemicalPotential"
                , chemicalPotential_function_value
                , "Return the chemical potential of the replica (if the moves\nsample a constant chemical potential ensemble)\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireMove::Replica::energy
        
            typedef ::SireUnits::Dimension::MolarEnergy ( ::SireMove::Replica::*energy_function_type)(  ) ;
            energy_function_type energy_function_value( &::SireMove::Replica::energy );
            
            Replica_exposer.def( 
                "energy"
                , energy_function_value
                , "Return the total energy of this replica" );
        
        }
        { //::SireMove::Replica::energyComponent
        
            typedef ::SireCAS::Symbol const & ( ::SireMove::Replica::*energyComponent_function_type)(  ) const;
            energyComponent_function_type energyComponent_function_value( &::SireMove::Replica::energyComponent );
            
            Replica_exposer.def( 
                "energyComponent"
                , energyComponent_function_value
                , bp::return_value_policy<bp::clone_const_reference>()
                , "Return the energy component that describes the Hamiltonian\nthat is sampled by this replica" );
        
        }
        { //::SireMove::Replica::ensemble
        
            typedef ::SireMove::Ensemble const & ( ::SireMove::Replica::*ensemble_function_type)(  ) const;
            ensemble_function_type ensemble_function_value( &::SireMove::Replica::ensemble );
            
            Replica_exposer.def( 
                "ensemble"
                , ensemble_function_value
                , bp::return_value_policy<bp::clone_const_reference>()
                , "Return the ensemble defined by the moves of this replica" );
        
        }
        { //::SireMove::Replica::fugacity
        
            typedef ::SireUnits::Dimension::Pressure ( ::SireMove::Replica::*fugacity_function_type)(  ) const;
            fugacity_function_type fugacity_function_value( &::SireMove::Replica::fugacity );
            
            Replica_exposer.def( 
                "fugacity"
                , fugacity_function_value
                , "Return the fugacity of the replica (if the moves sample\na constant chemical potential ensemble)\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireMove::Replica::isConstantChemicalPotential
        
            typedef bool ( ::SireMove::Replica::*isConstantChemicalPotential_function_type)(  ) const;
            isConstantChemicalPotential_function_type isConstantChemicalPotential_function_value( &::SireMove::Replica::isConstantChemicalPotential );
            
            Replica_exposer.def( 
                "isConstantChemicalPotential"
                , isConstantChemicalPotential_function_value
                , "Return whether the moves keep the same chemical potential" );
        
        }
        { //::SireMove::Replica::isConstantEnergy
        
            typedef bool ( ::SireMove::Replica::*isConstantEnergy_function_type)(  ) const;
            isConstantEnergy_function_type isConstantEnergy_function_value( &::SireMove::Replica::isConstantEnergy );
            
            Replica_exposer.def( 
                "isConstantEnergy"
                , isConstantEnergy_function_value
                , "Return whether or not this is a constant energy replica\n(all moves sample the same total energy)" );
        
        }
        { //::SireMove::Replica::isConstantFugacity
        
            typedef bool ( ::SireMove::Replica::*isConstantFugacity_function_type)(  ) const;
            isConstantFugacity_function_type isConstantFugacity_function_value( &::SireMove::Replica::isConstantFugacity );
            
            Replica_exposer.def( 
                "isConstantFugacity"
                , isConstantFugacity_function_value
                , "Return whether the moves keep the same fugacity" );
        
        }
        { //::SireMove::Replica::isConstantLambda
        
            typedef bool ( ::SireMove::Replica::*isConstantLambda_function_type)( ::SireCAS::Symbol const & ) const;
            isConstantLambda_function_type isConstantLambda_function_value( &::SireMove::Replica::isConstantLambda );
            
            Replica_exposer.def( 
                "isConstantLambda"
                , isConstantLambda_function_value
                , ( bp::arg("lam") )
                , "Return whether or not this is a constant lambda replica\n(there is a lambda component, lam, and all moves sample the\nsame value of this lambda coordinate)" );
        
        }
        { //::SireMove::Replica::isConstantNParticles
        
            typedef bool ( ::SireMove::Replica::*isConstantNParticles_function_type)(  ) const;
            isConstantNParticles_function_type isConstantNParticles_function_value( &::SireMove::Replica::isConstantNParticles );
            
            Replica_exposer.def( 
                "isConstantNParticles"
                , isConstantNParticles_function_value
                , "Return whether or not the number of particles is constant\n(all moves keep the same number of particles)" );
        
        }
        { //::SireMove::Replica::isConstantPressure
        
            typedef bool ( ::SireMove::Replica::*isConstantPressure_function_type)(  ) const;
            isConstantPressure_function_type isConstantPressure_function_value( &::SireMove::Replica::isConstantPressure );
            
            Replica_exposer.def( 
                "isConstantPressure"
                , isConstantPressure_function_value
                , "Return whether or not this is a constant pressure replica\n(all moves sample the same pressure)" );
        
        }
        { //::SireMove::Replica::isConstantTemperature
        
            typedef bool ( ::SireMove::Replica::*isConstantTemperature_function_type)(  ) const;
            isConstantTemperature_function_type isConstantTemperature_function_value( &::SireMove::Replica::isConstantTemperature );
            
            Replica_exposer.def( 
                "isConstantTemperature"
                , isConstantTemperature_function_value
                , "Return whether or not this is a constant temperature replica\n(all moves sample the same temperature)" );
        
        }
        { //::SireMove::Replica::isConstantVolume
        
            typedef bool ( ::SireMove::Replica::*isConstantVolume_function_type)(  ) const;
            isConstantVolume_function_type isConstantVolume_function_value( &::SireMove::Replica::isConstantVolume );
            
            Replica_exposer.def( 
                "isConstantVolume"
                , isConstantVolume_function_value
                , "Return whether or not this is a constant volume replica\n(all moves sample the same volume)" );
        
        }
        { //::SireMove::Replica::lambdaComponent
        
            typedef ::SireCAS::Symbol const & ( ::SireMove::Replica::*lambdaComponent_function_type)(  ) const;
            lambdaComponent_function_type lambdaComponent_function_value( &::SireMove::Replica::lambdaComponent );
            
            Replica_exposer.def( 
                "lambdaComponent"
                , lambdaComponent_function_value
                , bp::return_value_policy<bp::clone_const_reference>()
                , "Return the component that can be used to change the Hamiltonian\nfor Hamiltonian replica exchange - this is a null symbol if\nthis replica is not used in Hamiltonian replica exchange" );
        
        }
        { //::SireMove::Replica::lambdaValue
        
            typedef double ( ::SireMove::Replica::*lambdaValue_function_type)(  ) const;
            lambdaValue_function_type lambdaValue_function_value( &::SireMove::Replica::lambdaValue );
            
            Replica_exposer.def( 
                "lambdaValue"
                , lambdaValue_function_value
                , "Return the value of the lambda component - this returns 0 if\nthis replica is not suited for Hamiltonian replica exchange" );
        
        }
        Replica_exposer.def( bp::self != bp::self );
        { //::SireMove::Replica::operator=
        
            typedef ::SireMove::Replica & ( ::SireMove::Replica::*assign_function_type)( ::SireMove::Replica const & ) ;
            assign_function_type assign_function_value( &::SireMove::Replica::operator= );
            
            Replica_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        Replica_exposer.def( bp::self == bp::self );
        { //::SireMove::Replica::pressure
        
            typedef ::SireUnits::Dimension::Pressure ( ::SireMove::Replica::*pressure_function_type)(  ) const;
            pressure_function_type pressure_function_value( &::SireMove::Replica::pressure );
            
            Replica_exposer.def( 
                "pressure"
                , pressure_function_value
                , "Return the pressure of the replica (if the moves sample\na constant pressure ensemble)\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireMove::Replica::spaceProperty
        
            typedef ::SireBase::PropertyName const & ( ::SireMove::Replica::*spaceProperty_function_type)(  ) const;
            spaceProperty_function_type spaceProperty_function_value( &::SireMove::Replica::spaceProperty );
            
            Replica_exposer.def( 
                "spaceProperty"
                , spaceProperty_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the name of the property containing the simulation box\nthat is sampled by this replica - this is used to get the\nvolume of the simulation space" );
        
        }
        { //::SireMove::Replica::temperature
        
            typedef ::SireUnits::Dimension::Temperature ( ::SireMove::Replica::*temperature_function_type)(  ) const;
            temperature_function_type temperature_function_value( &::SireMove::Replica::temperature );
            
            Replica_exposer.def( 
                "temperature"
                , temperature_function_value
                , "Return the temperature of the replica (if the moves sample\na constant temperature ensemble)\nThrow: SireError::incompatible_error\n" );
        
        }
        { //::SireMove::Replica::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMove::Replica::typeName );
            
            Replica_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireMove::Replica::volume
        
            typedef ::SireUnits::Dimension::Volume ( ::SireMove::Replica::*volume_function_type)(  ) const;
            volume_function_type volume_function_value( &::SireMove::Replica::volume );
            
            Replica_exposer.def( 
                "volume"
                , volume_function_value
                , "Return the current volume of the replicas simulation space\n(this could be infinite)" );
        
        }
        Replica_exposer.staticmethod( "typeName" );
        Replica_exposer.def( "__copy__", &__copy__);
        Replica_exposer.def( "__deepcopy__", &__copy__);
        Replica_exposer.def( "clone", &__copy__);
        Replica_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMove::Replica >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Replica_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMove::Replica >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Replica_exposer.def_pickle(sire_pickle_suite< ::SireMove::Replica >());
        Replica_exposer.def( "__str__", &__str__< ::SireMove::Replica > );
        Replica_exposer.def( "__repr__", &__str__< ::SireMove::Replica > );
    }

}
