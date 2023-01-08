// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "IDOrSet_MGID_.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/incremint.h"

#include "SireError/errors.h"

#include "SireMol/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "mgid.h"

#include "mgidentifier.h"

#include "mgidx.h"

#include "mgname.h"

#include "mgnum.h"

#include "moleculegroups.h"

#include "tostring.h"

#include "mgid.h"

SireID::IDOrSet<SireMol::MGID> __copy__(const SireID::IDOrSet<SireMol::MGID> &other){ return SireID::IDOrSet<SireMol::MGID>(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_IDOrSet_MGID__class(){

    { //::SireID::IDOrSet< SireMol::MGID >
        typedef bp::class_< SireID::IDOrSet< SireMol::MGID >, bp::bases< SireMol::MGID, SireID::ID > > IDOrSet_MGID__exposer_t;
        IDOrSet_MGID__exposer_t IDOrSet_MGID__exposer = IDOrSet_MGID__exposer_t( "IDOrSet_MGID_", "", bp::init< >("") );
        bp::scope IDOrSet_MGID__scope( IDOrSet_MGID__exposer );
        IDOrSet_MGID__exposer.def( bp::init< SireMol::MGID const & >(( bp::arg("id") ), "") );
        IDOrSet_MGID__exposer.def( bp::init< SireMol::MGID const &, SireMol::MGID const & >(( bp::arg("id0"), bp::arg("id1") ), "") );
        IDOrSet_MGID__exposer.def( bp::init< QList< SireMol::MGIdentifier > const & >(( bp::arg("ids") ), "") );
        IDOrSet_MGID__exposer.def( bp::init< SireID::IDOrSet< SireMol::MGID > const & >(( bp::arg("ids") ), "") );
        IDOrSet_MGID__exposer.def( bp::init< SireID::IDOrSet< SireMol::MGID > const & >(( bp::arg("other") ), "") );
        { //::SireID::IDOrSet< SireMol::MGID >::IDs
        
            typedef SireID::IDOrSet< SireMol::MGID > exported_class_t;
            typedef ::QSet< SireMol::MGIdentifier > const & ( ::SireID::IDOrSet< SireMol::MGID >::*IDs_function_type)(  ) const;
            IDs_function_type IDs_function_value( &::SireID::IDOrSet< SireMol::MGID >::IDs );
            
            IDOrSet_MGID__exposer.def( 
                "IDs"
                , IDs_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireID::IDOrSet< SireMol::MGID >::hash
        
            typedef SireID::IDOrSet< SireMol::MGID > exported_class_t;
            typedef ::uint ( ::SireID::IDOrSet< SireMol::MGID >::*hash_function_type)(  ) const;
            hash_function_type hash_function_value( &::SireID::IDOrSet< SireMol::MGID >::hash );
            
            IDOrSet_MGID__exposer.def( 
                "hash"
                , hash_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireID::IDOrSet< SireMol::MGID >::isNull
        
            typedef SireID::IDOrSet< SireMol::MGID > exported_class_t;
            typedef bool ( ::SireID::IDOrSet< SireMol::MGID >::*isNull_function_type)(  ) const;
            isNull_function_type isNull_function_value( &::SireID::IDOrSet< SireMol::MGID >::isNull );
            
            IDOrSet_MGID__exposer.def( 
                "isNull"
                , isNull_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireID::IDOrSet< SireMol::MGID >::map
        
            typedef SireID::IDOrSet< SireMol::MGID > exported_class_t;
            typedef ::QList< SireMol::MGNum > ( ::SireID::IDOrSet< SireMol::MGID >::*map_function_type)( ::SireID::IDOrSet< SireMol::MGID >::SearchObject const & ) const;
            map_function_type map_function_value( &::SireID::IDOrSet< SireMol::MGID >::map );
            
            IDOrSet_MGID__exposer.def( 
                "map"
                , map_function_value
                , ( bp::arg("obj") )
                , bp::release_gil_policy()
                , "" );
        
        }
        IDOrSet_MGID__exposer.def( bp::self != bp::other< SireID::ID >() );
        IDOrSet_MGID__exposer.def( bp::self != bp::self );
        IDOrSet_MGID__exposer.def( bp::self != bp::other< SireMol::MGID >() );
        { //::SireID::IDOrSet< SireMol::MGID >::operator=
        
            typedef SireID::IDOrSet< SireMol::MGID > exported_class_t;
            typedef ::SireID::IDOrSet< SireMol::MGID > & ( ::SireID::IDOrSet< SireMol::MGID >::*assign_function_type)( ::SireID::IDOrSet< SireMol::MGID > const & ) ;
            assign_function_type assign_function_value( &::SireID::IDOrSet< SireMol::MGID >::operator= );
            
            IDOrSet_MGID__exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireID::IDOrSet< SireMol::MGID >::operator=
        
            typedef SireID::IDOrSet< SireMol::MGID > exported_class_t;
            typedef ::SireID::IDOrSet< SireMol::MGID > & ( ::SireID::IDOrSet< SireMol::MGID >::*assign_function_type)( ::SireMol::MGID const & ) ;
            assign_function_type assign_function_value( &::SireID::IDOrSet< SireMol::MGID >::operator= );
            
            IDOrSet_MGID__exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        IDOrSet_MGID__exposer.def( bp::self == bp::other< SireID::ID >() );
        IDOrSet_MGID__exposer.def( bp::self == bp::self );
        IDOrSet_MGID__exposer.def( bp::self == bp::other< SireMol::MGID >() );
        { //::SireID::IDOrSet< SireMol::MGID >::toString
        
            typedef SireID::IDOrSet< SireMol::MGID > exported_class_t;
            typedef ::QString ( ::SireID::IDOrSet< SireMol::MGID >::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireID::IDOrSet< SireMol::MGID >::toString );
            
            IDOrSet_MGID__exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireID::IDOrSet< SireMol::MGID >::typeName
        
            typedef SireID::IDOrSet< SireMol::MGID > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireID::IDOrSet< SireMol::MGID >::typeName );
            
            IDOrSet_MGID__exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireID::IDOrSet< SireMol::MGID >::what
        
            typedef SireID::IDOrSet< SireMol::MGID > exported_class_t;
            typedef char const * ( ::SireID::IDOrSet< SireMol::MGID >::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireID::IDOrSet< SireMol::MGID >::what );
            
            IDOrSet_MGID__exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        IDOrSet_MGID__exposer.staticmethod( "typeName" );
        IDOrSet_MGID__exposer.def( "__copy__", &__copy__);
        IDOrSet_MGID__exposer.def( "__deepcopy__", &__copy__);
        IDOrSet_MGID__exposer.def( "clone", &__copy__);
        IDOrSet_MGID__exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireID::IDOrSet<SireMol::MGID> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        IDOrSet_MGID__exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireID::IDOrSet<SireMol::MGID> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        IDOrSet_MGID__exposer.def_pickle(sire_pickle_suite< ::SireID::IDOrSet<SireMol::MGID> >());
        IDOrSet_MGID__exposer.def( "__str__", &__str__< ::SireID::IDOrSet<SireMol::MGID> > );
        IDOrSet_MGID__exposer.def( "__repr__", &__str__< ::SireID::IDOrSet<SireMol::MGID> > );
        IDOrSet_MGID__exposer.def( "__hash__", &::SireID::IDOrSet<SireMol::MGID>::hash );
    }

}
