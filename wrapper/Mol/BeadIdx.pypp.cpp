// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "BeadIdx.pypp.hpp"

namespace bp = boost::python;

#include "beadidx.h"

#include "beadidx.h"

#include "SireStream/datastream.h"

#include "beadid.h"

#include "beadidx.h"

#include "beadnum.h"

#include "beadid.h"

SireMol::BeadIdx __copy__(const SireMol::BeadIdx &other){ return SireMol::BeadIdx(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_BeadIdx_class(){

    { //::SireMol::BeadIdx
        typedef bp::class_< SireMol::BeadIdx, bp::bases< SireMol::BeadID, SireID::ID, SireID::IndexBase > > BeadIdx_exposer_t;
        BeadIdx_exposer_t BeadIdx_exposer = BeadIdx_exposer_t( "BeadIdx", "This is an ID object that is used to index CutGroups\n\nAuthor: Christopher Woods\n", bp::init< >("") );
        bp::scope BeadIdx_scope( BeadIdx_exposer );
        BeadIdx_exposer.def( bp::init< qint32 >(( bp::arg("idx") ), "") );
        BeadIdx_exposer.def( bp::init< SireMol::BeadIdx const & >(( bp::arg("other") ), "") );
        { //::SireMol::BeadIdx::hash
        
            typedef ::uint ( ::SireMol::BeadIdx::*hash_function_type)(  ) const;
            hash_function_type hash_function_value( &::SireMol::BeadIdx::hash );
            
            BeadIdx_exposer.def( 
                "hash"
                , hash_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::BeadIdx::isNull
        
            typedef bool ( ::SireMol::BeadIdx::*isNull_function_type)(  ) const;
            isNull_function_type isNull_function_value( &::SireMol::BeadIdx::isNull );
            
            BeadIdx_exposer.def( 
                "isNull"
                , isNull_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::BeadIdx::null
        
            typedef ::SireMol::BeadIdx ( *null_function_type )(  );
            null_function_type null_function_value( &::SireMol::BeadIdx::null );
            
            BeadIdx_exposer.def( 
                "null"
                , null_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        BeadIdx_exposer.def( bp::self != bp::other< SireID::ID >() );
        { //::SireMol::BeadIdx::operator=
        
            typedef ::SireMol::BeadIdx & ( ::SireMol::BeadIdx::*assign_function_type)( ::SireMol::BeadIdx const & ) ;
            assign_function_type assign_function_value( &::SireMol::BeadIdx::operator= );
            
            BeadIdx_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::BeadIdx::toString
        
            typedef ::QString ( ::SireMol::BeadIdx::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::BeadIdx::toString );
            
            BeadIdx_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::BeadIdx::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::BeadIdx::typeName );
            
            BeadIdx_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::BeadIdx::what
        
            typedef char const * ( ::SireMol::BeadIdx::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMol::BeadIdx::what );
            
            BeadIdx_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        BeadIdx_exposer.staticmethod( "null" );
        BeadIdx_exposer.staticmethod( "typeName" );
        BeadIdx_exposer.def( "__copy__", &__copy__);
        BeadIdx_exposer.def( "__deepcopy__", &__copy__);
        BeadIdx_exposer.def( "clone", &__copy__);
        BeadIdx_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::BeadIdx >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        BeadIdx_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::BeadIdx >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        BeadIdx_exposer.def_pickle(sire_pickle_suite< ::SireMol::BeadIdx >());
        BeadIdx_exposer.def( "__str__", &__str__< ::SireMol::BeadIdx > );
        BeadIdx_exposer.def( "__repr__", &__str__< ::SireMol::BeadIdx > );
        BeadIdx_exposer.def( "__hash__", &::SireMol::BeadIdx::hash );
    }

}
