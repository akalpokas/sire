// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "MonitorProperty.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireFF/ff.h"

#include "SireMol/errors.h"

#include "SireMol/moleculegroup.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "monitorproperty.h"

#include "system.h"

#include <QFile>

#include <QTextStream>

#include "monitorproperty.h"

SireSystem::MonitorProperty __copy__(const SireSystem::MonitorProperty &other){ return SireSystem::MonitorProperty(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_MonitorProperty_class(){

    { //::SireSystem::MonitorProperty
        typedef bp::class_< SireSystem::MonitorProperty, bp::bases< SireSystem::SystemMonitor, SireBase::Property > > MonitorProperty_exposer_t;
        MonitorProperty_exposer_t MonitorProperty_exposer = MonitorProperty_exposer_t( "MonitorProperty", bp::init< >() );
        bp::scope MonitorProperty_scope( MonitorProperty_exposer );
        MonitorProperty_exposer.def( bp::init< QString const & >(( bp::arg("property") )) );
        MonitorProperty_exposer.def( bp::init< QString const &, SireMol::MoleculeGroup const & >(( bp::arg("property"), bp::arg("molgroup") )) );
        MonitorProperty_exposer.def( bp::init< QString const &, SireMol::MGID const & >(( bp::arg("property"), bp::arg("mgid") )) );
        MonitorProperty_exposer.def( bp::init< QString const &, SireFF::FF const & >(( bp::arg("property"), bp::arg("forcefield") )) );
        MonitorProperty_exposer.def( bp::init< QString const &, SireFF::FFID const & >(( bp::arg("property"), bp::arg("ffid") )) );
        MonitorProperty_exposer.def( bp::init< SireSystem::MonitorProperty const & >(( bp::arg("other") )) );
        { //::SireSystem::MonitorProperty::clearStatistics
        
            typedef void ( ::SireSystem::MonitorProperty::*clearStatistics_function_type)(  ) ;
            clearStatistics_function_type clearStatistics_function_value( &::SireSystem::MonitorProperty::clearStatistics );
            
            MonitorProperty_exposer.def( 
                "clearStatistics"
                , clearStatistics_function_value );
        
        }
        { //::SireSystem::MonitorProperty::ffID
        
            typedef ::SireFF::FFID const & ( ::SireSystem::MonitorProperty::*ffID_function_type)(  ) const;
            ffID_function_type ffID_function_value( &::SireSystem::MonitorProperty::ffID );
            
            MonitorProperty_exposer.def( 
                "ffID"
                , ffID_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireSystem::MonitorProperty::mgID
        
            typedef ::SireMol::MGID const & ( ::SireSystem::MonitorProperty::*mgID_function_type)(  ) const;
            mgID_function_type mgID_function_value( &::SireSystem::MonitorProperty::mgID );
            
            MonitorProperty_exposer.def( 
                "mgID"
                , mgID_function_value
                , bp::return_value_policy<bp::clone_const_reference>() );
        
        }
        { //::SireSystem::MonitorProperty::monitor
        
            typedef void ( ::SireSystem::MonitorProperty::*monitor_function_type)( ::SireSystem::System & ) ;
            monitor_function_type monitor_function_value( &::SireSystem::MonitorProperty::monitor );
            
            MonitorProperty_exposer.def( 
                "monitor"
                , monitor_function_value
                , ( bp::arg("system") ) );
        
        }
        { //::SireSystem::MonitorProperty::monitoredMolecules
        
            typedef ::QList< SireMol::MolNum > ( ::SireSystem::MonitorProperty::*monitoredMolecules_function_type)(  ) const;
            monitoredMolecules_function_type monitoredMolecules_function_value( &::SireSystem::MonitorProperty::monitoredMolecules );
            
            MonitorProperty_exposer.def( 
                "monitoredMolecules"
                , monitoredMolecules_function_value );
        
        }
        { //::SireSystem::MonitorProperty::monitoringForceFieldProperty
        
            typedef bool ( ::SireSystem::MonitorProperty::*monitoringForceFieldProperty_function_type)(  ) const;
            monitoringForceFieldProperty_function_type monitoringForceFieldProperty_function_value( &::SireSystem::MonitorProperty::monitoringForceFieldProperty );
            
            MonitorProperty_exposer.def( 
                "monitoringForceFieldProperty"
                , monitoringForceFieldProperty_function_value );
        
        }
        { //::SireSystem::MonitorProperty::monitoringMoleculeProperty
        
            typedef bool ( ::SireSystem::MonitorProperty::*monitoringMoleculeProperty_function_type)(  ) const;
            monitoringMoleculeProperty_function_type monitoringMoleculeProperty_function_value( &::SireSystem::MonitorProperty::monitoringMoleculeProperty );
            
            MonitorProperty_exposer.def( 
                "monitoringMoleculeProperty"
                , monitoringMoleculeProperty_function_value );
        
        }
        { //::SireSystem::MonitorProperty::monitoringSystemProperty
        
            typedef bool ( ::SireSystem::MonitorProperty::*monitoringSystemProperty_function_type)(  ) const;
            monitoringSystemProperty_function_type monitoringSystemProperty_function_value( &::SireSystem::MonitorProperty::monitoringSystemProperty );
            
            MonitorProperty_exposer.def( 
                "monitoringSystemProperty"
                , monitoringSystemProperty_function_value );
        
        }
        MonitorProperty_exposer.def( bp::self != bp::self );
        { //::SireSystem::MonitorProperty::operator=
        
            typedef ::SireSystem::MonitorProperty & ( ::SireSystem::MonitorProperty::*assign_function_type)( ::SireSystem::MonitorProperty const & ) ;
            assign_function_type assign_function_value( &::SireSystem::MonitorProperty::operator= );
            
            MonitorProperty_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        MonitorProperty_exposer.def( bp::self == bp::self );
        { //::SireSystem::MonitorProperty::properties
        
            typedef ::QVector< SireBase::PropPtr< SireBase::Property > > ( ::SireSystem::MonitorProperty::*properties_function_type)(  ) const;
            properties_function_type properties_function_value( &::SireSystem::MonitorProperty::properties );
            
            MonitorProperty_exposer.def( 
                "properties"
                , properties_function_value );
        
        }
        { //::SireSystem::MonitorProperty::properties
        
            typedef ::QVector< SireBase::PropPtr< SireBase::Property > > ( ::SireSystem::MonitorProperty::*properties_function_type)( ::SireMol::MolNum ) const;
            properties_function_type properties_function_value( &::SireSystem::MonitorProperty::properties );
            
            MonitorProperty_exposer.def( 
                "properties"
                , properties_function_value
                , ( bp::arg("molnum") ) );
        
        }
        { //::SireSystem::MonitorProperty::property
        
            typedef ::QString const & ( ::SireSystem::MonitorProperty::*property_function_type)(  ) const;
            property_function_type property_function_value( &::SireSystem::MonitorProperty::property );
            
            MonitorProperty_exposer.def( 
                "property"
                , property_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireSystem::MonitorProperty::toString
        
            typedef ::QString ( ::SireSystem::MonitorProperty::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireSystem::MonitorProperty::toString );
            
            MonitorProperty_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireSystem::MonitorProperty::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireSystem::MonitorProperty::typeName );
            
            MonitorProperty_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireSystem::MonitorProperty::writeToDisk
        
            typedef void ( ::SireSystem::MonitorProperty::*writeToDisk_function_type)( ::QString const & ) ;
            writeToDisk_function_type writeToDisk_function_value( &::SireSystem::MonitorProperty::writeToDisk );
            
            MonitorProperty_exposer.def( 
                "writeToDisk"
                , writeToDisk_function_value
                , ( bp::arg("filename") ) );
        
        }
        MonitorProperty_exposer.staticmethod( "typeName" );
        MonitorProperty_exposer.def( "__copy__", &__copy__);
        MonitorProperty_exposer.def( "__deepcopy__", &__copy__);
        MonitorProperty_exposer.def( "clone", &__copy__);
        MonitorProperty_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireSystem::MonitorProperty >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        MonitorProperty_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireSystem::MonitorProperty >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        MonitorProperty_exposer.def( "__str__", &__str__< ::SireSystem::MonitorProperty > );
        MonitorProperty_exposer.def( "__repr__", &__str__< ::SireSystem::MonitorProperty > );
    }

}
