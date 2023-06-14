// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "CGAtomIdx.pypp.hpp"

namespace bp = boost::python;

#include "cgatomidx.h"

#include "molinfo.h"

#include "cgatomidx.h"

#include "cgidx.h"

#include "SireID/index.h"

SireMol::CGAtomIdx __copy__(const SireMol::CGAtomIdx &other){ return SireMol::CGAtomIdx(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_CGAtomIdx_class(){

    { //::SireMol::CGAtomIdx
        typedef bp::class_< SireMol::CGAtomIdx, bp::bases< SireMol::AtomID, SireID::ID > > CGAtomIdx_exposer_t;
        CGAtomIdx_exposer_t CGAtomIdx_exposer = CGAtomIdx_exposer_t( "CGAtomIdx", "This is the basic type used to ID atoms within a molecule. This\nprovides the fastest way of indexing atoms and is the base\ntype that all other AtomID classes map to.\n\nAuthor: Christopher Woods\n", bp::init< >("") );
        bp::scope CGAtomIdx_scope( CGAtomIdx_exposer );
        CGAtomIdx_exposer.def( bp::init< SireMol::CGIdx, SireID::Index >(( bp::arg("cgid"), bp::arg("atomid") ), "") );
        CGAtomIdx_exposer.def( bp::init< SireMol::CGAtomIdx const & >(( bp::arg("other") ), "") );
        { //::SireMol::CGAtomIdx::atom
        
            typedef ::SireID::Index ( ::SireMol::CGAtomIdx::*atom_function_type)(  ) const;
            atom_function_type atom_function_value( &::SireMol::CGAtomIdx::atom );
            
            CGAtomIdx_exposer.def( 
                "atom"
                , atom_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::CGAtomIdx::cutGroup
        
            typedef ::SireMol::CGIdx ( ::SireMol::CGAtomIdx::*cutGroup_function_type)(  ) const;
            cutGroup_function_type cutGroup_function_value( &::SireMol::CGAtomIdx::cutGroup );
            
            CGAtomIdx_exposer.def( 
                "cutGroup"
                , cutGroup_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::CGAtomIdx::hash
        
            typedef ::uint ( ::SireMol::CGAtomIdx::*hash_function_type)(  ) const;
            hash_function_type hash_function_value( &::SireMol::CGAtomIdx::hash );
            
            CGAtomIdx_exposer.def( 
                "hash"
                , hash_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::CGAtomIdx::isNull
        
            typedef bool ( ::SireMol::CGAtomIdx::*isNull_function_type)(  ) const;
            isNull_function_type isNull_function_value( &::SireMol::CGAtomIdx::isNull );
            
            CGAtomIdx_exposer.def( 
                "isNull"
                , isNull_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::CGAtomIdx::map
        
            typedef ::QList< SireMol::AtomIdx > ( ::SireMol::CGAtomIdx::*map_function_type)( ::SireMol::MolInfo const & ) const;
            map_function_type map_function_value( &::SireMol::CGAtomIdx::map );
            
            CGAtomIdx_exposer.def( 
                "map"
                , map_function_value
                , ( bp::arg("molinfo") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::CGAtomIdx::null
        
            typedef ::SireMol::CGAtomIdx ( *null_function_type )(  );
            null_function_type null_function_value( &::SireMol::CGAtomIdx::null );
            
            CGAtomIdx_exposer.def( 
                "null"
                , null_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        CGAtomIdx_exposer.def( bp::self != bp::other< SireID::ID >() );
        CGAtomIdx_exposer.def( bp::self != bp::self );
        { //::SireMol::CGAtomIdx::operator=
        
            typedef ::SireMol::CGAtomIdx & ( ::SireMol::CGAtomIdx::*assign_function_type)( ::SireMol::CGAtomIdx const & ) ;
            assign_function_type assign_function_value( &::SireMol::CGAtomIdx::operator= );
            
            CGAtomIdx_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        CGAtomIdx_exposer.def( bp::self == bp::other< SireID::ID >() );
        CGAtomIdx_exposer.def( bp::self == bp::self );
        { //::SireMol::CGAtomIdx::toString
        
            typedef ::QString ( ::SireMol::CGAtomIdx::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::CGAtomIdx::toString );
            
            CGAtomIdx_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::CGAtomIdx::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::CGAtomIdx::typeName );
            
            CGAtomIdx_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::CGAtomIdx::what
        
            typedef char const * ( ::SireMol::CGAtomIdx::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMol::CGAtomIdx::what );
            
            CGAtomIdx_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        CGAtomIdx_exposer.staticmethod( "null" );
        CGAtomIdx_exposer.staticmethod( "typeName" );
        CGAtomIdx_exposer.def( "__copy__", &__copy__);
        CGAtomIdx_exposer.def( "__deepcopy__", &__copy__);
        CGAtomIdx_exposer.def( "clone", &__copy__);
        CGAtomIdx_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::CGAtomIdx >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CGAtomIdx_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::CGAtomIdx >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        CGAtomIdx_exposer.def_pickle(sire_pickle_suite< ::SireMol::CGAtomIdx >());
        CGAtomIdx_exposer.def( "__str__", &__str__< ::SireMol::CGAtomIdx > );
        CGAtomIdx_exposer.def( "__repr__", &__str__< ::SireMol::CGAtomIdx > );
        CGAtomIdx_exposer.def( "__hash__", &::SireMol::CGAtomIdx::hash );
    }

}
