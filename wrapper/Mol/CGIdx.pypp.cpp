// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "CGIdx.pypp.hpp"

namespace bp = boost::python;

#include "cgidx.h"

#include "cgidx.h"

#include "SireID/index.h"

#include "cgatomidx.h"

SireMol::CGIdx __copy__(const SireMol::CGIdx &other){ return SireMol::CGIdx(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_CGIdx_class(){

    { //::SireMol::CGIdx
        typedef bp::class_< SireMol::CGIdx, bp::bases< SireMol::CGID, SireID::ID, SireID::IndexBase > > CGIdx_exposer_t;
        CGIdx_exposer_t CGIdx_exposer = CGIdx_exposer_t( "CGIdx", "This is an ID object that is used to index CutGroups\n\nAuthor: Christopher Woods\n", bp::init< >("") );
        bp::scope CGIdx_scope( CGIdx_exposer );
        CGIdx_exposer.def( bp::init< qint32 >(( bp::arg("idx") ), "") );
        CGIdx_exposer.def( bp::init< SireMol::CGIdx const & >(( bp::arg("other") ), "") );
        { //::SireMol::CGIdx::hash
        
            typedef ::uint ( ::SireMol::CGIdx::*hash_function_type)(  ) const;
            hash_function_type hash_function_value( &::SireMol::CGIdx::hash );
            
            CGIdx_exposer.def( 
                "hash"
                , hash_function_value
                , "" );
        
        }
        { //::SireMol::CGIdx::isNull
        
            typedef bool ( ::SireMol::CGIdx::*isNull_function_type)(  ) const;
            isNull_function_type isNull_function_value( &::SireMol::CGIdx::isNull );
            
            CGIdx_exposer.def( 
                "isNull"
                , isNull_function_value
                , "" );
        
        }
        { //::SireMol::CGIdx::map
        
            typedef ::QList< SireMol::CGIdx > ( ::SireMol::CGIdx::*map_function_type)( ::SireMol::MolInfo const & ) const;
            map_function_type map_function_value( &::SireMol::CGIdx::map );
            
            CGIdx_exposer.def( 
                "map"
                , map_function_value
                , ( bp::arg("molinfo") )
                , "" );
        
        }
        { //::SireMol::CGIdx::null
        
            typedef ::SireMol::CGIdx ( *null_function_type )(  );
            null_function_type null_function_value( &::SireMol::CGIdx::null );
            
            CGIdx_exposer.def( 
                "null"
                , null_function_value
                , "" );
        
        }
        CGIdx_exposer.def( bp::self + bp::other< SireID::Index >() );
        { //::SireMol::CGIdx::operator=
        
            typedef ::SireMol::CGIdx & ( ::SireMol::CGIdx::*assign_function_type)( ::SireMol::CGIdx const & ) ;
            assign_function_type assign_function_value( &::SireMol::CGIdx::operator= );
            
            CGIdx_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::CGIdx::toString
        
            typedef ::QString ( ::SireMol::CGIdx::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::CGIdx::toString );
            
            CGIdx_exposer.def( 
                "toString"
                , toString_function_value
                , "" );
        
        }
        { //::SireMol::CGIdx::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::CGIdx::typeName );
            
            CGIdx_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireMol::CGIdx::what
        
            typedef char const * ( ::SireMol::CGIdx::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMol::CGIdx::what );
            
            CGIdx_exposer.def( 
                "what"
                , what_function_value
                , "" );
        
        }
        CGIdx_exposer.staticmethod( "null" );
        CGIdx_exposer.staticmethod( "typeName" );
        CGIdx_exposer.def( other<SireID::Index>() + self );
        CGIdx_exposer.def( "__copy__", &__copy__);
        CGIdx_exposer.def( "__deepcopy__", &__copy__);
        CGIdx_exposer.def( "clone", &__copy__);
        CGIdx_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::CGIdx >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CGIdx_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::CGIdx >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CGIdx_exposer.def( "__setstate__", &__setstate__base64< ::SireMol::CGIdx > );
        CGIdx_exposer.def( "__getstate__", &__getstate__base64< ::SireMol::CGIdx > );
        CGIdx_exposer.def( "__str__", &__str__< ::SireMol::CGIdx > );
        CGIdx_exposer.def( "__repr__", &__str__< ::SireMol::CGIdx > );
        CGIdx_exposer.def( "__hash__", &::SireMol::CGIdx::hash );
    }

}
