// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "ID.pypp.hpp"

namespace bp = boost::python;

#include "id.h"

#include "name.h"

#include "number.h"

#include <QDebug>

#include "id.h"

#include "Helpers/str.hpp"

void register_ID_class(){

    { //::SireID::ID
        typedef bp::class_< SireID::ID, boost::noncopyable > ID_exposer_t;
        ID_exposer_t ID_exposer = ID_exposer_t( "ID", bp::no_init );
        bp::scope ID_scope( ID_exposer );
        { //::SireID::ID::hash
        
            typedef ::uint ( ::SireID::ID::*hash_function_type)(  ) const;
            hash_function_type hash_function_value( &::SireID::ID::hash );
            
            ID_exposer.def( 
                "hash"
                , hash_function_value );
        
        }
        { //::SireID::ID::isNull
        
            typedef bool ( ::SireID::ID::*isNull_function_type)(  ) const;
            isNull_function_type isNull_function_value( &::SireID::ID::isNull );
            
            ID_exposer.def( 
                "isNull"
                , isNull_function_value );
        
        }
        ID_exposer.def( bp::self != bp::self );
        ID_exposer.def( bp::self == bp::self );
        { //::SireID::ID::toString
        
            typedef ::QString ( ::SireID::ID::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireID::ID::toString );
            
            ID_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireID::ID::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireID::ID::typeName );
            
            ID_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireID::ID::what
        
            typedef char const * ( ::SireID::ID::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireID::ID::what );
            
            ID_exposer.def( 
                "what"
                , what_function_value );
        
        }
        ID_exposer.staticmethod( "typeName" );
        ID_exposer.def( "__str__", &__str__< ::SireID::ID > );
        ID_exposer.def( "__repr__", &__str__< ::SireID::ID > );
        ID_exposer.def( "__hash__", &::SireID::ID::hash );
    }

}
