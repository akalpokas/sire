// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "IntegerArrayProperty.pypp.hpp"

namespace bp = boost::python;

#include "SireError/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "arrayproperty.hpp"

#include "booleanproperty.h"

#include "numberproperty.h"

#include "propertylist.h"

#include "stringproperty.h"

#include "tostring.h"

#include "variantproperty.h"

#include "propertylist.h"

SireBase::IntegerArrayProperty __copy__(const SireBase::IntegerArrayProperty &other){ return SireBase::IntegerArrayProperty(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_IntegerArrayProperty_class(){

    { //::SireBase::IntegerArrayProperty
        typedef bp::class_< SireBase::IntegerArrayProperty, bp::bases< SireBase::ArrayProperty<long long>, SireBase::Property > > IntegerArrayProperty_exposer_t;
        IntegerArrayProperty_exposer_t IntegerArrayProperty_exposer = IntegerArrayProperty_exposer_t( "IntegerArrayProperty", "", bp::init< >("") );
        bp::scope IntegerArrayProperty_scope( IntegerArrayProperty_exposer );
        IntegerArrayProperty_exposer.def( bp::init< QList< long long > const & >(( bp::arg("array") ), "") );
        IntegerArrayProperty_exposer.def( bp::init< QVector< long long > const & >(( bp::arg("array") ), "") );
        IntegerArrayProperty_exposer.def( bp::init< SireBase::StringArrayProperty const & >(( bp::arg("array") ), "") );
        IntegerArrayProperty_exposer.def( bp::init< SireBase::DoubleArrayProperty const & >(( bp::arg("array") ), "") );
        IntegerArrayProperty_exposer.def( bp::init< SireBase::PropertyList const & >(( bp::arg("array") ), "") );
        IntegerArrayProperty_exposer.def( bp::init< SireBase::IntegerArrayProperty const & >(( bp::arg("other") ), "") );
        { //::SireBase::IntegerArrayProperty::asABoolean
        
            typedef bool ( ::SireBase::IntegerArrayProperty::*asABoolean_function_type)(  ) const;
            asABoolean_function_type asABoolean_function_value( &::SireBase::IntegerArrayProperty::asABoolean );
            
            IntegerArrayProperty_exposer.def( 
                "asABoolean"
                , asABoolean_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireBase::IntegerArrayProperty::asADouble
        
            typedef double ( ::SireBase::IntegerArrayProperty::*asADouble_function_type)(  ) const;
            asADouble_function_type asADouble_function_value( &::SireBase::IntegerArrayProperty::asADouble );
            
            IntegerArrayProperty_exposer.def( 
                "asADouble"
                , asADouble_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireBase::IntegerArrayProperty::asAString
        
            typedef ::QString ( ::SireBase::IntegerArrayProperty::*asAString_function_type)(  ) const;
            asAString_function_type asAString_function_value( &::SireBase::IntegerArrayProperty::asAString );
            
            IntegerArrayProperty_exposer.def( 
                "asAString"
                , asAString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireBase::IntegerArrayProperty::asAnArray
        
            typedef ::SireBase::PropertyList ( ::SireBase::IntegerArrayProperty::*asAnArray_function_type)(  ) const;
            asAnArray_function_type asAnArray_function_value( &::SireBase::IntegerArrayProperty::asAnArray );
            
            IntegerArrayProperty_exposer.def( 
                "asAnArray"
                , asAnArray_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireBase::IntegerArrayProperty::asAnInteger
        
            typedef int ( ::SireBase::IntegerArrayProperty::*asAnInteger_function_type)(  ) const;
            asAnInteger_function_type asAnInteger_function_value( &::SireBase::IntegerArrayProperty::asAnInteger );
            
            IntegerArrayProperty_exposer.def( 
                "asAnInteger"
                , asAnInteger_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireBase::IntegerArrayProperty::isABoolean
        
            typedef bool ( ::SireBase::IntegerArrayProperty::*isABoolean_function_type)(  ) const;
            isABoolean_function_type isABoolean_function_value( &::SireBase::IntegerArrayProperty::isABoolean );
            
            IntegerArrayProperty_exposer.def( 
                "isABoolean"
                , isABoolean_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireBase::IntegerArrayProperty::isADouble
        
            typedef bool ( ::SireBase::IntegerArrayProperty::*isADouble_function_type)(  ) const;
            isADouble_function_type isADouble_function_value( &::SireBase::IntegerArrayProperty::isADouble );
            
            IntegerArrayProperty_exposer.def( 
                "isADouble"
                , isADouble_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireBase::IntegerArrayProperty::isAString
        
            typedef bool ( ::SireBase::IntegerArrayProperty::*isAString_function_type)(  ) const;
            isAString_function_type isAString_function_value( &::SireBase::IntegerArrayProperty::isAString );
            
            IntegerArrayProperty_exposer.def( 
                "isAString"
                , isAString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireBase::IntegerArrayProperty::isAnInteger
        
            typedef bool ( ::SireBase::IntegerArrayProperty::*isAnInteger_function_type)(  ) const;
            isAnInteger_function_type isAnInteger_function_value( &::SireBase::IntegerArrayProperty::isAnInteger );
            
            IntegerArrayProperty_exposer.def( 
                "isAnInteger"
                , isAnInteger_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        IntegerArrayProperty_exposer.def( bp::self != bp::self );
        IntegerArrayProperty_exposer.def( bp::self + bp::self );
        { //::SireBase::IntegerArrayProperty::operator=
        
            typedef ::SireBase::IntegerArrayProperty & ( ::SireBase::IntegerArrayProperty::*assign_function_type)( ::SireBase::IntegerArrayProperty const & ) ;
            assign_function_type assign_function_value( &::SireBase::IntegerArrayProperty::operator= );
            
            IntegerArrayProperty_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        IntegerArrayProperty_exposer.def( bp::self == bp::self );
        { //::SireBase::IntegerArrayProperty::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireBase::IntegerArrayProperty::typeName );
            
            IntegerArrayProperty_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        IntegerArrayProperty_exposer.staticmethod( "typeName" );
        IntegerArrayProperty_exposer.def( "__copy__", &__copy__);
        IntegerArrayProperty_exposer.def( "__deepcopy__", &__copy__);
        IntegerArrayProperty_exposer.def( "clone", &__copy__);
        IntegerArrayProperty_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireBase::IntegerArrayProperty >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        IntegerArrayProperty_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireBase::IntegerArrayProperty >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        IntegerArrayProperty_exposer.def_pickle(sire_pickle_suite< ::SireBase::IntegerArrayProperty >());
        IntegerArrayProperty_exposer.def( "__str__", &__str__< ::SireBase::IntegerArrayProperty > );
        IntegerArrayProperty_exposer.def( "__repr__", &__str__< ::SireBase::IntegerArrayProperty > );
        IntegerArrayProperty_exposer.def( "__len__", &__len_size< ::SireBase::IntegerArrayProperty > );
        IntegerArrayProperty_exposer.def( "__getitem__", &::SireBase::IntegerArrayProperty::getitem );
    }

}
