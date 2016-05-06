// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "AtomIdx.pypp.hpp"

namespace bp = boost::python;

#include "atomidx.h"

#include "atomidx.h"

#include "SireMol/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "atom.h"

#include "atomid.h"

#include "atomidentifier.h"

#include "chain.h"

#include "cutgroup.h"

#include "editor.hpp"

#include "groupatomids.h"

#include "molatomid.h"

#include "moleculegroup.h"

#include "moleculegroups.h"

#include "molecules.h"

#include "molinfo.h"

#include "mover.hpp"

#include "residue.h"

#include "segment.h"

#include "selector.hpp"

#include "tostring.h"

#include "withatoms.h"

#include <QDebug>

#include "atomid.h"

SireMol::AtomIdx __copy__(const SireMol::AtomIdx &other){ return SireMol::AtomIdx(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_AtomIdx_class(){

    { //::SireMol::AtomIdx
        typedef bp::class_< SireMol::AtomIdx, bp::bases< SireMol::AtomID, SireID::ID, SireID::IndexBase > > AtomIdx_exposer_t;
        AtomIdx_exposer_t AtomIdx_exposer = AtomIdx_exposer_t( "AtomIdx", bp::init< >() );
        bp::scope AtomIdx_scope( AtomIdx_exposer );
        AtomIdx_exposer.def( bp::init< qint32 >(( bp::arg("idx") )) );
        AtomIdx_exposer.def( bp::init< SireMol::AtomIdx const & >(( bp::arg("other") )) );
        { //::SireMol::AtomIdx::hash
        
            typedef ::uint ( ::SireMol::AtomIdx::*hash_function_type)(  ) const;
            hash_function_type hash_function_value( &::SireMol::AtomIdx::hash );
            
            AtomIdx_exposer.def( 
                "hash"
                , hash_function_value );
        
        }
        { //::SireMol::AtomIdx::isNull
        
            typedef bool ( ::SireMol::AtomIdx::*isNull_function_type)(  ) const;
            isNull_function_type isNull_function_value( &::SireMol::AtomIdx::isNull );
            
            AtomIdx_exposer.def( 
                "isNull"
                , isNull_function_value );
        
        }
        { //::SireMol::AtomIdx::map
        
            typedef ::QList< SireMol::AtomIdx > ( ::SireMol::AtomIdx::*map_function_type)( ::SireMol::MolInfo const & ) const;
            map_function_type map_function_value( &::SireMol::AtomIdx::map );
            
            AtomIdx_exposer.def( 
                "map"
                , map_function_value
                , ( bp::arg("molinfo") ) );
        
        }
        { //::SireMol::AtomIdx::null
        
            typedef ::SireMol::AtomIdx ( *null_function_type )(  );
            null_function_type null_function_value( &::SireMol::AtomIdx::null );
            
            AtomIdx_exposer.def( 
                "null"
                , null_function_value );
        
        }
        { //::SireMol::AtomIdx::operator=
        
            typedef ::SireMol::AtomIdx & ( ::SireMol::AtomIdx::*assign_function_type)( ::SireMol::AtomIdx const & ) ;
            assign_function_type assign_function_value( &::SireMol::AtomIdx::operator= );
            
            AtomIdx_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        { //::SireMol::AtomIdx::toString
        
            typedef ::QString ( ::SireMol::AtomIdx::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::AtomIdx::toString );
            
            AtomIdx_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireMol::AtomIdx::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::AtomIdx::typeName );
            
            AtomIdx_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireMol::AtomIdx::what
        
            typedef char const * ( ::SireMol::AtomIdx::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMol::AtomIdx::what );
            
            AtomIdx_exposer.def( 
                "what"
                , what_function_value );
        
        }
        AtomIdx_exposer.staticmethod( "null" );
        AtomIdx_exposer.staticmethod( "typeName" );
        AtomIdx_exposer.def( "__copy__", &__copy__);
        AtomIdx_exposer.def( "__deepcopy__", &__copy__);
        AtomIdx_exposer.def( "clone", &__copy__);
        AtomIdx_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::AtomIdx >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AtomIdx_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::AtomIdx >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AtomIdx_exposer.def( "__str__", &__str__< ::SireMol::AtomIdx > );
        AtomIdx_exposer.def( "__repr__", &__str__< ::SireMol::AtomIdx > );
        AtomIdx_exposer.def( "__hash__", &::SireMol::AtomIdx::hash );
    }

}
