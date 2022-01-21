// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "SegChainID.pypp.hpp"

namespace bp = boost::python;

#include "SireStream/datastream.h"

#include "groupgroupids.h"

#include "groupgroupids.h"

SireMol::GroupGroupID<SireMol::SegID, SireMol::ChainID> __copy__(const SireMol::GroupGroupID<SireMol::SegID, SireMol::ChainID> &other){ return SireMol::GroupGroupID<SireMol::SegID, SireMol::ChainID>(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_SegChainID_class(){

    { //::SireMol::GroupGroupID< SireMol::SegID, SireMol::ChainID >
        typedef bp::class_< SireMol::GroupGroupID< SireMol::SegID, SireMol::ChainID >, bp::bases< SireMol::GroupAtomIDBase, SireMol::AtomID, SireID::ID > > SegChainID_exposer_t;
        SegChainID_exposer_t SegChainID_exposer = SegChainID_exposer_t( "SegChainID", "", bp::init< >("") );
        bp::scope SegChainID_scope( SegChainID_exposer );
        SegChainID_exposer.def( bp::init< SireMol::SegID const &, SireMol::ChainID const & >(( bp::arg("group0"), bp::arg("group1") ), "") );
        SegChainID_exposer.def( bp::init< SireMol::GroupGroupID< SireMol::SegID, SireMol::ChainID > const & >(( bp::arg("other") ), "") );
        { //::SireMol::GroupGroupID< SireMol::SegID, SireMol::ChainID >::hash
        
            typedef SireMol::GroupGroupID< SireMol::SegID, SireMol::ChainID > exported_class_t;
            typedef ::uint ( ::SireMol::GroupGroupID< SireMol::SegID, SireMol::ChainID >::*hash_function_type)(  ) const;
            hash_function_type hash_function_value( &::SireMol::GroupGroupID< SireMol::SegID, SireMol::ChainID >::hash );
            
            SegChainID_exposer.def( 
                "hash"
                , hash_function_value
                , "" );
        
        }
        { //::SireMol::GroupGroupID< SireMol::SegID, SireMol::ChainID >::isNull
        
            typedef SireMol::GroupGroupID< SireMol::SegID, SireMol::ChainID > exported_class_t;
            typedef bool ( ::SireMol::GroupGroupID< SireMol::SegID, SireMol::ChainID >::*isNull_function_type)(  ) const;
            isNull_function_type isNull_function_value( &::SireMol::GroupGroupID< SireMol::SegID, SireMol::ChainID >::isNull );
            
            SegChainID_exposer.def( 
                "isNull"
                , isNull_function_value
                , "" );
        
        }
        { //::SireMol::GroupGroupID< SireMol::SegID, SireMol::ChainID >::map
        
            typedef SireMol::GroupGroupID< SireMol::SegID, SireMol::ChainID > exported_class_t;
            typedef ::QList< SireMol::AtomIdx > ( ::SireMol::GroupGroupID< SireMol::SegID, SireMol::ChainID >::*map_function_type)( ::SireMol::MolInfo const & ) const;
            map_function_type map_function_value( &::SireMol::GroupGroupID< SireMol::SegID, SireMol::ChainID >::map );
            
            SegChainID_exposer.def( 
                "map"
                , map_function_value
                , ( bp::arg("molinfo") )
                , "" );
        
        }
        SegChainID_exposer.def( bp::self != bp::self );
        SegChainID_exposer.def( bp::self == bp::other< SireID::ID >() );
        SegChainID_exposer.def( bp::self == bp::self );
        { //::SireMol::GroupGroupID< SireMol::SegID, SireMol::ChainID >::toString
        
            typedef SireMol::GroupGroupID< SireMol::SegID, SireMol::ChainID > exported_class_t;
            typedef ::QString ( ::SireMol::GroupGroupID< SireMol::SegID, SireMol::ChainID >::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::GroupGroupID< SireMol::SegID, SireMol::ChainID >::toString );
            
            SegChainID_exposer.def( 
                "toString"
                , toString_function_value
                , "" );
        
        }
        { //::SireMol::GroupGroupID< SireMol::SegID, SireMol::ChainID >::typeName
        
            typedef SireMol::GroupGroupID< SireMol::SegID, SireMol::ChainID > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::GroupGroupID< SireMol::SegID, SireMol::ChainID >::typeName );
            
            SegChainID_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireMol::GroupGroupID< SireMol::SegID, SireMol::ChainID >::what
        
            typedef SireMol::GroupGroupID< SireMol::SegID, SireMol::ChainID > exported_class_t;
            typedef char const * ( ::SireMol::GroupGroupID< SireMol::SegID, SireMol::ChainID >::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMol::GroupGroupID< SireMol::SegID, SireMol::ChainID >::what );
            
            SegChainID_exposer.def( 
                "what"
                , what_function_value
                , "" );
        
        }
        SegChainID_exposer.staticmethod( "typeName" );
        SegChainID_exposer.def( "__copy__", &__copy__);
        SegChainID_exposer.def( "__deepcopy__", &__copy__);
        SegChainID_exposer.def( "clone", &__copy__);
        SegChainID_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::GroupGroupID<SireMol::SegID, SireMol::ChainID> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        SegChainID_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::GroupGroupID<SireMol::SegID, SireMol::ChainID> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        SegChainID_exposer.def( "__setstate__", &__setstate__base64< ::SireMol::GroupGroupID<SireMol::SegID, SireMol::ChainID> > );
        SegChainID_exposer.def( "__getstate__", &__getstate__base64< ::SireMol::GroupGroupID<SireMol::SegID, SireMol::ChainID> > );
        SegChainID_exposer.def( "__str__", &__str__< ::SireMol::GroupGroupID<SireMol::SegID, SireMol::ChainID> > );
        SegChainID_exposer.def( "__repr__", &__str__< ::SireMol::GroupGroupID<SireMol::SegID, SireMol::ChainID> > );
        SegChainID_exposer.def( "__hash__", &::SireMol::GroupGroupID<SireMol::SegID, SireMol::ChainID>::hash );
    }

}
