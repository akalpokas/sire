// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "Name.pypp.hpp"

namespace bp = boost::python;

#include "name.h"

#include "name.h"

#include "Qt/qdatastream.hpp"

const char* pvt_get_name(const SireID::Name&){ return "SireID::Name";}

#include "Helpers/release_gil_policy.hpp"

void register_Name_class(){

    { //::SireID::Name
        typedef bp::class_< SireID::Name, boost::noncopyable > Name_exposer_t;
        Name_exposer_t Name_exposer = Name_exposer_t( "Name", "This is the base class of all Name ID objects. A Name is used\nto provide an object with a human-readable name that can be\nused to identify the object, e.g. identifying atoms by their\nname within a residue, or identifying forcefields by their name.\n\nA name does not have to uniquely identify an object, though it\nhelps, as generally only the first object with a specified name\nis returned if there are in fact multiple objects with the same\nname.\n\nAuthor: Christopher Woods\n", bp::no_init );
        bp::scope Name_scope( Name_exposer );
        { //::SireID::Name::hash
        
            typedef ::uint ( ::SireID::Name::*hash_function_type)(  ) const;
            hash_function_type hash_function_value( &::SireID::Name::hash );
            
            Name_exposer.def( 
                "hash"
                , hash_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireID::Name::isCaseSensitive
        
            typedef bool ( ::SireID::Name::*isCaseSensitive_function_type)(  ) const;
            isCaseSensitive_function_type isCaseSensitive_function_value( &::SireID::Name::isCaseSensitive );
            
            Name_exposer.def( 
                "isCaseSensitive"
                , isCaseSensitive_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireID::Name::isEmpty
        
            typedef bool ( ::SireID::Name::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireID::Name::isEmpty );
            
            Name_exposer.def( 
                "isEmpty"
                , isEmpty_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireID::Name::isNull
        
            typedef bool ( ::SireID::Name::*isNull_function_type)(  ) const;
            isNull_function_type isNull_function_value( &::SireID::Name::isNull );
            
            Name_exposer.def( 
                "isNull"
                , isNull_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireID::Name::value
        
            typedef ::QString const & ( ::SireID::Name::*value_function_type)(  ) const;
            value_function_type value_function_value( &::SireID::Name::value );
            
            Name_exposer.def( 
                "value"
                , value_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        Name_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireID::Name >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Name_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireID::Name >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Name_exposer.def_pickle(sire_pickle_suite< ::SireID::Name >());
        Name_exposer.def( "__str__", &pvt_get_name);
        Name_exposer.def( "__repr__", &pvt_get_name);
        Name_exposer.def( "__hash__", &::SireID::Name::hash );
    }

}
