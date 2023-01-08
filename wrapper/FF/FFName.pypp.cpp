// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "FFName.pypp.hpp"

namespace bp = boost::python;

#include "ffname.h"

#include "ffname.h"

#include "forcefields.h"

#include "ffidx.h"

SireFF::FFName __copy__(const SireFF::FFName &other){ return SireFF::FFName(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_FFName_class(){

    { //::SireFF::FFName
        typedef bp::class_< SireFF::FFName, bp::bases< SireFF::FFID, SireID::ID, SireID::Name > > FFName_exposer_t;
        FFName_exposer_t FFName_exposer = FFName_exposer_t( "FFName", "This class holds the name of a forcefield\n\nAuthor: Christopher Woods\n", bp::init< >("") );
        bp::scope FFName_scope( FFName_exposer );
        FFName_exposer.def( bp::init< QString const & >(( bp::arg("name") ), "") );
        FFName_exposer.def( bp::init< QString const &, SireID::CaseSensitivity >(( bp::arg("name"), bp::arg("case_sensitivity") ), "") );
        FFName_exposer.def( bp::init< SireFF::FFName const & >(( bp::arg("other") ), "") );
        { //::SireFF::FFName::hash
        
            typedef ::uint ( ::SireFF::FFName::*hash_function_type)(  ) const;
            hash_function_type hash_function_value( &::SireFF::FFName::hash );
            
            FFName_exposer.def( 
                "hash"
                , hash_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::FFName::isNull
        
            typedef bool ( ::SireFF::FFName::*isNull_function_type)(  ) const;
            isNull_function_type isNull_function_value( &::SireFF::FFName::isNull );
            
            FFName_exposer.def( 
                "isNull"
                , isNull_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::FFName::map
        
            typedef ::QList< SireFF::FFIdx > ( ::SireFF::FFName::*map_function_type)( ::SireFF::ForceFields const & ) const;
            map_function_type map_function_value( &::SireFF::FFName::map );
            
            FFName_exposer.def( 
                "map"
                , map_function_value
                , ( bp::arg("ffields") )
                , bp::release_gil_policy()
                , "Map this name to the index of the matching forcefield in the\npassed ForceFields object\nThrow: SireFF::missing_forcefield\n" );
        
        }
        FFName_exposer.def( bp::self != bp::self );
        { //::SireFF::FFName::operator=
        
            typedef ::SireFF::FFName & ( ::SireFF::FFName::*assign_function_type)( ::SireFF::FFName const & ) ;
            assign_function_type assign_function_value( &::SireFF::FFName::operator= );
            
            FFName_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        FFName_exposer.def( bp::self == bp::other< SireID::ID >() );
        FFName_exposer.def( bp::self == bp::self );
        { //::SireFF::FFName::toString
        
            typedef ::QString ( ::SireFF::FFName::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireFF::FFName::toString );
            
            FFName_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::FFName::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireFF::FFName::typeName );
            
            FFName_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireFF::FFName::what
        
            typedef char const * ( ::SireFF::FFName::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireFF::FFName::what );
            
            FFName_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        FFName_exposer.staticmethod( "typeName" );
        FFName_exposer.def( "__copy__", &__copy__);
        FFName_exposer.def( "__deepcopy__", &__copy__);
        FFName_exposer.def( "clone", &__copy__);
        FFName_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireFF::FFName >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        FFName_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireFF::FFName >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        FFName_exposer.def_pickle(sire_pickle_suite< ::SireFF::FFName >());
        FFName_exposer.def( "__str__", &__str__< ::SireFF::FFName > );
        FFName_exposer.def( "__repr__", &__str__< ::SireFF::FFName > );
        FFName_exposer.def( "__hash__", &::SireFF::FFName::hash );
    }

}
