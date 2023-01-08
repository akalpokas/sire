// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "MolIdx.pypp.hpp"

namespace bp = boost::python;

#include "molidx.h"

#include "molidx.h"

#include "molid.h"

#include "molidx.h"

#include "molnum.h"

#include "molname.h"

#include "moleculegroup.h"

#include "moleculegroups.h"

#include "mover.hpp"

SireMol::MolIdx __copy__(const SireMol::MolIdx &other){ return SireMol::MolIdx(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_MolIdx_class(){

    { //::SireMol::MolIdx
        typedef bp::class_< SireMol::MolIdx, bp::bases< SireMol::MolID, SireID::ID, SireID::IndexBase > > MolIdx_exposer_t;
        MolIdx_exposer_t MolIdx_exposer = MolIdx_exposer_t( "MolIdx", "This is an ID object that is used to index molecules (e.g. index\nin a list or array, or in a MoleculeGroup).\n\nAuthor: Christopher Woods\n", bp::init< >("") );
        bp::scope MolIdx_scope( MolIdx_exposer );
        MolIdx_exposer.def( bp::init< qint32 >(( bp::arg("idx") ), "") );
        MolIdx_exposer.def( bp::init< SireMol::MolIdx const & >(( bp::arg("other") ), "") );
        { //::SireMol::MolIdx::hash
        
            typedef ::uint ( ::SireMol::MolIdx::*hash_function_type)(  ) const;
            hash_function_type hash_function_value( &::SireMol::MolIdx::hash );
            
            MolIdx_exposer.def( 
                "hash"
                , hash_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::MolIdx::isNull
        
            typedef bool ( ::SireMol::MolIdx::*isNull_function_type)(  ) const;
            isNull_function_type isNull_function_value( &::SireMol::MolIdx::isNull );
            
            MolIdx_exposer.def( 
                "isNull"
                , isNull_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::MolIdx::map
        
            typedef ::QList< SireMol::MolNum > ( ::SireMol::MolIdx::*map_function_type)( ::SireMol::Molecules const & ) const;
            map_function_type map_function_value( &::SireMol::MolIdx::map );
            
            MolIdx_exposer.def( 
                "map"
                , map_function_value
                , ( bp::arg("molecules") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::MolIdx::map
        
            typedef ::QList< SireMol::MolNum > ( ::SireMol::MolIdx::*map_function_type)( ::SireMol::MoleculeGroup const & ) const;
            map_function_type map_function_value( &::SireMol::MolIdx::map );
            
            MolIdx_exposer.def( 
                "map"
                , map_function_value
                , ( bp::arg("molgroup") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::MolIdx::map
        
            typedef ::QList< SireMol::MolNum > ( ::SireMol::MolIdx::*map_function_type)( ::SireMol::MolGroupsBase const & ) const;
            map_function_type map_function_value( &::SireMol::MolIdx::map );
            
            MolIdx_exposer.def( 
                "map"
                , map_function_value
                , ( bp::arg("molgroups") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::MolIdx::null
        
            typedef ::SireMol::MolIdx ( *null_function_type )(  );
            null_function_type null_function_value( &::SireMol::MolIdx::null );
            
            MolIdx_exposer.def( 
                "null"
                , null_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::MolIdx::operator=
        
            typedef ::SireMol::MolIdx & ( ::SireMol::MolIdx::*assign_function_type)( ::SireMol::MolIdx const & ) ;
            assign_function_type assign_function_value( &::SireMol::MolIdx::operator= );
            
            MolIdx_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::MolIdx::toString
        
            typedef ::QString ( ::SireMol::MolIdx::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::MolIdx::toString );
            
            MolIdx_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::MolIdx::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::MolIdx::typeName );
            
            MolIdx_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::MolIdx::what
        
            typedef char const * ( ::SireMol::MolIdx::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMol::MolIdx::what );
            
            MolIdx_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        MolIdx_exposer.staticmethod( "null" );
        MolIdx_exposer.staticmethod( "typeName" );
        MolIdx_exposer.def( "__copy__", &__copy__);
        MolIdx_exposer.def( "__deepcopy__", &__copy__);
        MolIdx_exposer.def( "clone", &__copy__);
        MolIdx_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::MolIdx >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        MolIdx_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::MolIdx >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        MolIdx_exposer.def_pickle(sire_pickle_suite< ::SireMol::MolIdx >());
        MolIdx_exposer.def( "__str__", &__str__< ::SireMol::MolIdx > );
        MolIdx_exposer.def( "__repr__", &__str__< ::SireMol::MolIdx > );
        MolIdx_exposer.def( "__hash__", &::SireMol::MolIdx::hash );
    }

}
