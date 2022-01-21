// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "ResProp.pypp.hpp"

namespace bp = boost::python;

#include "resproperty.hpp"

#include "resproperty.hpp"

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_ResProp_class(){

    { //::SireMol::ResProp
        typedef bp::class_< SireMol::ResProp, bp::bases< SireMol::MolViewProperty, SireBase::Property >, boost::noncopyable > ResProp_exposer_t;
        ResProp_exposer_t ResProp_exposer = ResProp_exposer_t( "ResProp", "Small class used to provide a common base for all ResProperty types", bp::no_init );
        bp::scope ResProp_scope( ResProp_exposer );
        { //::SireMol::ResProp::assertCanConvert
        
            typedef void ( ::SireMol::ResProp::*assertCanConvert_function_type)( ::QVariant const & ) const;
            assertCanConvert_function_type assertCanConvert_function_value( &::SireMol::ResProp::assertCanConvert );
            
            ResProp_exposer.def( 
                "assertCanConvert"
                , assertCanConvert_function_value
                , ( bp::arg("value") )
                , "" );
        
        }
        { //::SireMol::ResProp::assignFrom
        
            typedef void ( ::SireMol::ResProp::*assignFrom_function_type)( ::SireMol::ResProperty< QVariant > const & ) ;
            assignFrom_function_type assignFrom_function_value( &::SireMol::ResProp::assignFrom );
            
            ResProp_exposer.def( 
                "assignFrom"
                , assignFrom_function_value
                , ( bp::arg("values") )
                , "" );
        
        }
        { //::SireMol::ResProp::canConvert
        
            typedef bool ( ::SireMol::ResProp::*canConvert_function_type)( ::QVariant const & ) const;
            canConvert_function_type canConvert_function_value( &::SireMol::ResProp::canConvert );
            
            ResProp_exposer.def( 
                "canConvert"
                , canConvert_function_value
                , ( bp::arg("value") )
                , "" );
        
        }
        { //::SireMol::ResProp::toVariant
        
            typedef ::SireMol::ResProperty< QVariant > ( ::SireMol::ResProp::*toVariant_function_type)(  ) const;
            toVariant_function_type toVariant_function_value( &::SireMol::ResProp::toVariant );
            
            ResProp_exposer.def( 
                "toVariant"
                , toVariant_function_value
                , "" );
        
        }
        ResProp_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::ResProp >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ResProp_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::ResProp >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ResProp_exposer.def( "__setstate__", &__setstate__base64< ::SireMol::ResProp > );
        ResProp_exposer.def( "__getstate__", &__getstate__base64< ::SireMol::ResProp > );
        ResProp_exposer.def( "__str__", &__str__< ::SireMol::ResProp > );
        ResProp_exposer.def( "__repr__", &__str__< ::SireMol::ResProp > );
    }

}
