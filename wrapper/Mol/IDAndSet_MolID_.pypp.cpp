// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "IDAndSet_MolID_.pypp.hpp"

namespace bp = boost::python;

#include "SireBase/incremint.h"

#include "SireError/errors.h"

#include "SireMol/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "atomidx.h"

#include "mgidx.h"

#include "molatomid.h"

#include "moleculegroup.h"

#include "moleculegroups.h"

#include "molecules.h"

#include "molid.h"

#include "molidx.h"

#include "molname.h"

#include "molnum.h"

#include "mover.hpp"

#include "specifymol.h"

#include "tostring.h"

#include "molid.h"

SireID::IDAndSet<SireMol::MolID> __copy__(const SireID::IDAndSet<SireMol::MolID> &other){ return SireID::IDAndSet<SireMol::MolID>(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_IDAndSet_MolID__class(){

    { //::SireID::IDAndSet< SireMol::MolID >
        typedef bp::class_< SireID::IDAndSet< SireMol::MolID >, bp::bases< SireMol::MolID, SireID::ID > > IDAndSet_MolID__exposer_t;
        IDAndSet_MolID__exposer_t IDAndSet_MolID__exposer = IDAndSet_MolID__exposer_t( "IDAndSet_MolID_", bp::init< >() );
        bp::scope IDAndSet_MolID__scope( IDAndSet_MolID__exposer );
        IDAndSet_MolID__exposer.def( bp::init< SireMol::MolID const & >(( bp::arg("id") )) );
        IDAndSet_MolID__exposer.def( bp::init< SireMol::MolID const &, SireMol::MolID const & >(( bp::arg("id0"), bp::arg("id1") )) );
        IDAndSet_MolID__exposer.def( bp::init< QList< SireMol::MolIdentifier > const & >(( bp::arg("ids") )) );
        IDAndSet_MolID__exposer.def( bp::init< SireID::IDAndSet< SireMol::MolID > const & >(( bp::arg("other") )) );
        { //::SireID::IDAndSet< SireMol::MolID >::IDs
        
            typedef SireID::IDAndSet< SireMol::MolID > exported_class_t;
            typedef ::QSet< SireMol::MolIdentifier > const & ( ::SireID::IDAndSet< SireMol::MolID >::*IDs_function_type)(  ) const;
            IDs_function_type IDs_function_value( &::SireID::IDAndSet< SireMol::MolID >::IDs );
            
            IDAndSet_MolID__exposer.def( 
                "IDs"
                , IDs_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireID::IDAndSet< SireMol::MolID >::hash
        
            typedef SireID::IDAndSet< SireMol::MolID > exported_class_t;
            typedef ::uint ( ::SireID::IDAndSet< SireMol::MolID >::*hash_function_type)(  ) const;
            hash_function_type hash_function_value( &::SireID::IDAndSet< SireMol::MolID >::hash );
            
            IDAndSet_MolID__exposer.def( 
                "hash"
                , hash_function_value );
        
        }
        { //::SireID::IDAndSet< SireMol::MolID >::isNull
        
            typedef SireID::IDAndSet< SireMol::MolID > exported_class_t;
            typedef bool ( ::SireID::IDAndSet< SireMol::MolID >::*isNull_function_type)(  ) const;
            isNull_function_type isNull_function_value( &::SireID::IDAndSet< SireMol::MolID >::isNull );
            
            IDAndSet_MolID__exposer.def( 
                "isNull"
                , isNull_function_value );
        
        }
        { //::SireID::IDAndSet< SireMol::MolID >::map
        
            typedef SireID::IDAndSet< SireMol::MolID > exported_class_t;
            typedef ::QList< SireMol::MolNum > ( ::SireID::IDAndSet< SireMol::MolID >::*map_function_type)( ::SireMol::Molecules const & ) const;
            map_function_type map_function_value( &::SireID::IDAndSet< SireMol::MolID >::map );
            
            IDAndSet_MolID__exposer.def( 
                "map"
                , map_function_value
                , ( bp::arg("mols") ) );
        
        }
        { //::SireID::IDAndSet< SireMol::MolID >::map
        
            typedef SireID::IDAndSet< SireMol::MolID > exported_class_t;
            typedef ::QList< SireMol::MolNum > ( ::SireID::IDAndSet< SireMol::MolID >::*map_function_type)( ::SireMol::MoleculeGroup const & ) const;
            map_function_type map_function_value( &::SireID::IDAndSet< SireMol::MolID >::map );
            
            IDAndSet_MolID__exposer.def( 
                "map"
                , map_function_value
                , ( bp::arg("molgroup") ) );
        
        }
        { //::SireID::IDAndSet< SireMol::MolID >::map
        
            typedef SireID::IDAndSet< SireMol::MolID > exported_class_t;
            typedef ::QList< SireMol::MolNum > ( ::SireID::IDAndSet< SireMol::MolID >::*map_function_type)( ::SireMol::MolGroupsBase const & ) const;
            map_function_type map_function_value( &::SireID::IDAndSet< SireMol::MolID >::map );
            
            IDAndSet_MolID__exposer.def( 
                "map"
                , map_function_value
                , ( bp::arg("molgroups") ) );
        
        }
        IDAndSet_MolID__exposer.def( bp::self != bp::other< SireID::ID >() );
        IDAndSet_MolID__exposer.def( bp::self != bp::self );
        IDAndSet_MolID__exposer.def( bp::self != bp::other< SireMol::MolID >() );
        { //::SireID::IDAndSet< SireMol::MolID >::operator=
        
            typedef SireID::IDAndSet< SireMol::MolID > exported_class_t;
            typedef ::SireID::IDAndSet< SireMol::MolID > & ( ::SireID::IDAndSet< SireMol::MolID >::*assign_function_type)( ::SireID::IDAndSet< SireMol::MolID > const & ) ;
            assign_function_type assign_function_value( &::SireID::IDAndSet< SireMol::MolID >::operator= );
            
            IDAndSet_MolID__exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        { //::SireID::IDAndSet< SireMol::MolID >::operator=
        
            typedef SireID::IDAndSet< SireMol::MolID > exported_class_t;
            typedef ::SireID::IDAndSet< SireMol::MolID > & ( ::SireID::IDAndSet< SireMol::MolID >::*assign_function_type)( ::SireMol::MolID const & ) ;
            assign_function_type assign_function_value( &::SireID::IDAndSet< SireMol::MolID >::operator= );
            
            IDAndSet_MolID__exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        IDAndSet_MolID__exposer.def( bp::self == bp::other< SireID::ID >() );
        IDAndSet_MolID__exposer.def( bp::self == bp::self );
        IDAndSet_MolID__exposer.def( bp::self == bp::other< SireMol::MolID >() );
        { //::SireID::IDAndSet< SireMol::MolID >::toString
        
            typedef SireID::IDAndSet< SireMol::MolID > exported_class_t;
            typedef ::QString ( ::SireID::IDAndSet< SireMol::MolID >::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireID::IDAndSet< SireMol::MolID >::toString );
            
            IDAndSet_MolID__exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireID::IDAndSet< SireMol::MolID >::typeName
        
            typedef SireID::IDAndSet< SireMol::MolID > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireID::IDAndSet< SireMol::MolID >::typeName );
            
            IDAndSet_MolID__exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        { //::SireID::IDAndSet< SireMol::MolID >::what
        
            typedef SireID::IDAndSet< SireMol::MolID > exported_class_t;
            typedef char const * ( ::SireID::IDAndSet< SireMol::MolID >::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireID::IDAndSet< SireMol::MolID >::what );
            
            IDAndSet_MolID__exposer.def( 
                "what"
                , what_function_value );
        
        }
        IDAndSet_MolID__exposer.staticmethod( "typeName" );
        IDAndSet_MolID__exposer.def( "__copy__", &__copy__);
        IDAndSet_MolID__exposer.def( "__deepcopy__", &__copy__);
        IDAndSet_MolID__exposer.def( "clone", &__copy__);
        IDAndSet_MolID__exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireID::IDAndSet<SireMol::MolID> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        IDAndSet_MolID__exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireID::IDAndSet<SireMol::MolID> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        IDAndSet_MolID__exposer.def( "__str__", &__str__< ::SireID::IDAndSet<SireMol::MolID> > );
        IDAndSet_MolID__exposer.def( "__repr__", &__str__< ::SireID::IDAndSet<SireMol::MolID> > );
        IDAndSet_MolID__exposer.def( "__hash__", &::SireID::IDAndSet<SireMol::MolID>::hash );
    }

}
