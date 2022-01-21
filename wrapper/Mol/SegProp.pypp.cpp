// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "SegProp.pypp.hpp"

namespace bp = boost::python;

#include "segproperty.hpp"

#include "segproperty.hpp"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_SegProp_class(){

    { //::SireMol::SegProp
        typedef bp::class_< SireMol::SegProp, bp::bases< SireMol::MolViewProperty, SireBase::Property >, boost::noncopyable > SegProp_exposer_t;
        SegProp_exposer_t SegProp_exposer = SegProp_exposer_t( "SegProp", "Small class used to provide a common base for all SegProperty types", bp::no_init );
        bp::scope SegProp_scope( SegProp_exposer );
        { //::SireMol::SegProp::assertCanConvert
        
            typedef void ( ::SireMol::SegProp::*assertCanConvert_function_type)( ::QVariant const & ) const;
            assertCanConvert_function_type assertCanConvert_function_value( &::SireMol::SegProp::assertCanConvert );
            
            SegProp_exposer.def( 
                "assertCanConvert"
                , assertCanConvert_function_value
                , ( bp::arg("value") )
                , "" );
        
        }
        { //::SireMol::SegProp::assignFrom
        
            typedef void ( ::SireMol::SegProp::*assignFrom_function_type)( ::SireMol::SegProperty< QVariant > const & ) ;
            assignFrom_function_type assignFrom_function_value( &::SireMol::SegProp::assignFrom );
            
            SegProp_exposer.def( 
                "assignFrom"
                , assignFrom_function_value
                , ( bp::arg("values") )
                , "" );
        
        }
        { //::SireMol::SegProp::canConvert
        
            typedef bool ( ::SireMol::SegProp::*canConvert_function_type)( ::QVariant const & ) const;
            canConvert_function_type canConvert_function_value( &::SireMol::SegProp::canConvert );
            
            SegProp_exposer.def( 
                "canConvert"
                , canConvert_function_value
                , ( bp::arg("value") )
                , "" );
        
        }
        { //::SireMol::SegProp::toVariant
        
            typedef ::SireMol::SegProperty< QVariant > ( ::SireMol::SegProp::*toVariant_function_type)(  ) const;
            toVariant_function_type toVariant_function_value( &::SireMol::SegProp::toVariant );
            
            SegProp_exposer.def( 
                "toVariant"
                , toVariant_function_value
                , "" );
        
        }
        SegProp_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::SegProp >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        SegProp_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::SegProp >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        SegProp_exposer.def( "__setstate__", &__setstate__base64< ::SireMol::SegProp > );
        SegProp_exposer.def( "__getstate__", &__getstate__base64< ::SireMol::SegProp > );
        SegProp_exposer.def( "__str__", &__str__< ::SireMol::SegProp > );
        SegProp_exposer.def( "__repr__", &__str__< ::SireMol::SegProp > );
    }

}
