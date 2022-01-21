// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Mopac.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/findexe.h"

#include "SireBase/sire_process.h"

#include "SireBase/tempdir.h"

#include "SireError/errors.h"

#include "SireMol/atomcharges.h"

#include "SireMol/atomcoords.h"

#include "SireMol/atomelements.h"

#include "SireMol/element.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "SireUnits/dimensions.h"

#include "SireUnits/units.h"

#include "mopac.h"

#include "qmpotential.h"

#include "tostring.h"

#include <QDebug>

#include "mopac.h"

Squire::Mopac __copy__(const Squire::Mopac &other){ return Squire::Mopac(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_Mopac_class(){

    { //::Squire::Mopac
        typedef bp::class_< Squire::Mopac, bp::bases< Squire::QMProgram, SireBase::Property > > Mopac_exposer_t;
        Mopac_exposer_t Mopac_exposer = Mopac_exposer_t( "Mopac", "This is a wrapper that allows Mopac to be used to calculate\nsemiempirical QM energies\n\nAuthor: Christopher Woods\n", bp::init< >("Constructor") );
        bp::scope Mopac_scope( Mopac_exposer );
        Mopac_exposer.def( bp::init< QString const & >(( bp::arg("mopac") ), "Construct, passing in the location of the Mopac executable") );
        Mopac_exposer.def( bp::init< Squire::Mopac const & >(( bp::arg("other") ), "Copy constructor") );
        { //::Squire::Mopac::calculateCharges
        
            typedef ::SireMol::AtomCharges ( ::Squire::Mopac::*calculateCharges_function_type)( ::SireMol::Molecule const &,::SireBase::PropertyMap const & ) const;
            calculateCharges_function_type calculateCharges_function_value( &::Squire::Mopac::calculateCharges );
            
            Mopac_exposer.def( 
                "calculateCharges"
                , calculateCharges_function_value
                , ( bp::arg("molecule"), bp::arg("map") )
                , "Use Mopac to calculate the partial charges for the passed molecule\n(using the total charge set in this program)" );
        
        }
        { //::Squire::Mopac::chargeTemplate
        
            typedef ::QString const & ( ::Squire::Mopac::*chargeTemplate_function_type)(  ) const;
            chargeTemplate_function_type chargeTemplate_function_value( &::Squire::Mopac::chargeTemplate );
            
            Mopac_exposer.def( 
                "chargeTemplate"
                , chargeTemplate_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the template for the command file to be used to get Mopac\nto calculate atomic partial charges." );
        
        }
        { //::Squire::Mopac::energyTemplate
        
            typedef ::QString const & ( ::Squire::Mopac::*energyTemplate_function_type)(  ) const;
            energyTemplate_function_type energyTemplate_function_value( &::Squire::Mopac::energyTemplate );
            
            Mopac_exposer.def( 
                "energyTemplate"
                , energyTemplate_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the template for the command file to be used to get Mopac\nto calculate the energy." );
        
        }
        { //::Squire::Mopac::environment
        
            typedef ::QHash< QString, QString > const & ( ::Squire::Mopac::*environment_function_type)(  ) const;
            environment_function_type environment_function_value( &::Squire::Mopac::environment );
            
            Mopac_exposer.def( 
                "environment"
                , environment_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return all of the environmental variables that are to be set explicitly\nwhen Mopac is run. This does not include any environmental variables\nthat have not been explicitly set, but do have values" );
        
        }
        { //::Squire::Mopac::environment
        
            typedef ::QString ( ::Squire::Mopac::*environment_function_type)( ::QString const & ) const;
            environment_function_type environment_function_value( &::Squire::Mopac::environment );
            
            Mopac_exposer.def( 
                "environment"
                , environment_function_value
                , ( bp::arg("variable") )
                , "Return the value of the explicitly set environmental variable variable.\nA null string is returned if this variable has not been set\nexplicitly (this does not mean the variable doesnt exist - merely\nthat a specific value has not been set)" );
        
        }
        { //::Squire::Mopac::executable
        
            typedef ::QString ( ::Squire::Mopac::*executable_function_type)(  ) const;
            executable_function_type executable_function_value( &::Squire::Mopac::executable );
            
            Mopac_exposer.def( 
                "executable"
                , executable_function_value
                , "Return the executable (full path and also arguments) to be used. This\nis null if the executable is searched for in the path" );
        
        }
        { //::Squire::Mopac::forceTemplate
        
            typedef ::QString const & ( ::Squire::Mopac::*forceTemplate_function_type)(  ) const;
            forceTemplate_function_type forceTemplate_function_value( &::Squire::Mopac::forceTemplate );
            
            Mopac_exposer.def( 
                "forceTemplate"
                , forceTemplate_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the template for the command file to be used to get Mopac\nto calculate the forces." );
        
        }
        { //::Squire::Mopac::method
        
            typedef ::QString const & ( ::Squire::Mopac::*method_function_type)(  ) const;
            method_function_type method_function_value( &::Squire::Mopac::method );
            
            Mopac_exposer.def( 
                "method"
                , method_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the QM method to be used to calculate the energy or\nforce (e.g. AM1, PM3). This will substitute for\n@QM_METHOD@ in the energy and force command file templates" );
        
        }
        { //::Squire::Mopac::mopacInputFilename
        
            typedef ::QString const & ( ::Squire::Mopac::*mopacInputFilename_function_type)(  ) const;
            mopacInputFilename_function_type mopacInputFilename_function_value( &::Squire::Mopac::mopacInputFilename );
            
            Mopac_exposer.def( 
                "mopacInputFilename"
                , mopacInputFilename_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the name of the file mopac is using for input. This\nis hardcoded into Mopac and depends on the system on which this\nprogram is running, e.g. FOR005 is the name for the Mopac\nI am using on OS X" );
        
        }
        { //::Squire::Mopac::mopacOutputFilename
        
            typedef ::QString const & ( ::Squire::Mopac::*mopacOutputFilename_function_type)(  ) const;
            mopacOutputFilename_function_type mopacOutputFilename_function_value( &::Squire::Mopac::mopacOutputFilename );
            
            Mopac_exposer.def( 
                "mopacOutputFilename"
                , mopacOutputFilename_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the name of the file mopac is using for output. This\nis hardcoded into Mopac and depends on the system on which this\nprogram is running, e.g. FOR006 is the name for the Mopac\nI am using on OS X" );
        
        }
        Mopac_exposer.def( bp::self != bp::self );
        { //::Squire::Mopac::operator=
        
            typedef ::Squire::Mopac & ( ::Squire::Mopac::*assign_function_type)( ::Squire::Mopac const & ) ;
            assign_function_type assign_function_value( &::Squire::Mopac::operator= );
            
            Mopac_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        Mopac_exposer.def( bp::self == bp::self );
        { //::Squire::Mopac::setChargeTemplate
        
            typedef void ( ::Squire::Mopac::*setChargeTemplate_function_type)( ::QString const & ) ;
            setChargeTemplate_function_type setChargeTemplate_function_value( &::Squire::Mopac::setChargeTemplate );
            
            Mopac_exposer.def( 
                "setChargeTemplate"
                , setChargeTemplate_function_value
                , ( bp::arg("charge_template") )
                , "Set the template for the command file to be used to get\nMopac to calculate atomic partial charges. The following tags will\n" );
        
        }
        { //::Squire::Mopac::setEnergyTemplate
        
            typedef void ( ::Squire::Mopac::*setEnergyTemplate_function_type)( ::QString const & ) ;
            setEnergyTemplate_function_type setEnergyTemplate_function_value( &::Squire::Mopac::setEnergyTemplate );
            
            Mopac_exposer.def( 
                "setEnergyTemplate"
                , setEnergyTemplate_function_value
                , ( bp::arg("energy_template") )
                , "Set the template for the command file to be used to get\nMopac to calculate an energy. The following tags will\n" );
        
        }
        { //::Squire::Mopac::setEnvironment
        
            typedef void ( ::Squire::Mopac::*setEnvironment_function_type)( ::QString const &,::QString const & ) ;
            setEnvironment_function_type setEnvironment_function_value( &::Squire::Mopac::setEnvironment );
            
            Mopac_exposer.def( 
                "setEnvironment"
                , setEnvironment_function_value
                , ( bp::arg("variable"), bp::arg("value") )
                , "Set the environmental variable variable to have the value value\nwhen the Mopac executable is run. This replaces any existing\nvalue of this environmental variable" );
        
        }
        { //::Squire::Mopac::setExecutable
        
            typedef void ( ::Squire::Mopac::*setExecutable_function_type)( ::QString const & ) ;
            setExecutable_function_type setExecutable_function_value( &::Squire::Mopac::setExecutable );
            
            Mopac_exposer.def( 
                "setExecutable"
                , setExecutable_function_value
                , ( bp::arg("mopac_exe") )
                , "Set the Mopac executable (full path and also arguments) to be used" );
        
        }
        { //::Squire::Mopac::setForceTemplate
        
            typedef void ( ::Squire::Mopac::*setForceTemplate_function_type)( ::QString const & ) ;
            setForceTemplate_function_type setForceTemplate_function_value( &::Squire::Mopac::setForceTemplate );
            
            Mopac_exposer.def( 
                "setForceTemplate"
                , setForceTemplate_function_value
                , ( bp::arg("force_template") )
                , "Set the template for the command file to be used to get\nMopac to calculate the forces. The following tags will\n" );
        
        }
        { //::Squire::Mopac::setMethod
        
            typedef void ( ::Squire::Mopac::*setMethod_function_type)( ::QString const & ) ;
            setMethod_function_type setMethod_function_value( &::Squire::Mopac::setMethod );
            
            Mopac_exposer.def( 
                "setMethod"
                , setMethod_function_value
                , ( bp::arg("method") )
                , "Set the QM method to be used to calculate the energy or\nforce (e.g. AM1, PM3). This will substitute for\n@QM_METHOD@ in the energy and force command file templates" );
        
        }
        { //::Squire::Mopac::setMopacInputFilename
        
            typedef void ( ::Squire::Mopac::*setMopacInputFilename_function_type)( ::QString const & ) ;
            setMopacInputFilename_function_type setMopacInputFilename_function_value( &::Squire::Mopac::setMopacInputFilename );
            
            Mopac_exposer.def( 
                "setMopacInputFilename"
                , setMopacInputFilename_function_value
                , ( bp::arg("filename") )
                , "Tell this interface which name Mopac uses for the input file. This\nis hardcoded into Mopac and depends on the system on which this\nprogram is running, e.g. FOR005 is the name for the Mopac\nI am using on OS X" );
        
        }
        { //::Squire::Mopac::setMopacOutputFilename
        
            typedef void ( ::Squire::Mopac::*setMopacOutputFilename_function_type)( ::QString const & ) ;
            setMopacOutputFilename_function_type setMopacOutputFilename_function_value( &::Squire::Mopac::setMopacOutputFilename );
            
            Mopac_exposer.def( 
                "setMopacOutputFilename"
                , setMopacOutputFilename_function_value
                , ( bp::arg("filename") )
                , "Tell this interface which name Mopac uses for the output file. This\nis hardcoded into Mopac and depends on the system on which this\nprogram is running, e.g. FOR006 is the name for the Mopac\nI am using on OS X" );
        
        }
        { //::Squire::Mopac::setTotalCharge
        
            typedef void ( ::Squire::Mopac::*setTotalCharge_function_type)( int ) ;
            setTotalCharge_function_type setTotalCharge_function_value( &::Squire::Mopac::setTotalCharge );
            
            Mopac_exposer.def( 
                "setTotalCharge"
                , setTotalCharge_function_value
                , ( bp::arg("charge") )
                , "Set the total charge of the system (in unit charges)" );
        
        }
        { //::Squire::Mopac::totalCharge
        
            typedef int ( ::Squire::Mopac::*totalCharge_function_type)(  ) const;
            totalCharge_function_type totalCharge_function_value( &::Squire::Mopac::totalCharge );
            
            Mopac_exposer.def( 
                "totalCharge"
                , totalCharge_function_value
                , "Return the total charge of the system" );
        
        }
        { //::Squire::Mopac::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::Squire::Mopac::typeName );
            
            Mopac_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        Mopac_exposer.staticmethod( "typeName" );
        Mopac_exposer.def( "__copy__", &__copy__);
        Mopac_exposer.def( "__deepcopy__", &__copy__);
        Mopac_exposer.def( "clone", &__copy__);
        Mopac_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::Squire::Mopac >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Mopac_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::Squire::Mopac >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Mopac_exposer.def( "__setstate__", &__setstate__base64< ::Squire::Mopac > );
        Mopac_exposer.def( "__getstate__", &__getstate__base64< ::Squire::Mopac > );
        Mopac_exposer.def( "__str__", &__str__< ::Squire::Mopac > );
        Mopac_exposer.def( "__repr__", &__str__< ::Squire::Mopac > );
    }

}
