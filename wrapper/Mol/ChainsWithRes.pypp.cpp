// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "ChainsWithRes.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/errors.h"

#include "SireStream/datastream.h"

#include "withres.h"

#include "withres.h"

SireMol::ChainsWithRes __copy__(const SireMol::ChainsWithRes &other){ return SireMol::ChainsWithRes(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

void register_ChainsWithRes_class(){

    { //::SireMol::ChainsWithRes
        typedef bp::class_< SireMol::ChainsWithRes, bp::bases< SireMol::ChainID, SireID::ID > > ChainsWithRes_exposer_t;
        ChainsWithRes_exposer_t ChainsWithRes_exposer = ChainsWithRes_exposer_t( "ChainsWithRes", "This ID class identifies chains that contain residues that\nmatch the passed ResID\n\nAuthor: Christopher Woods\n", bp::init< >("Null constructor") );
        bp::scope ChainsWithRes_scope( ChainsWithRes_exposer );
        ChainsWithRes_exposer.def( bp::init< SireMol::ResID const & >(( bp::arg("resid") ), "Construct from the passed residue ID") );
        ChainsWithRes_exposer.def( bp::init< SireMol::ChainsWithRes const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMol::ChainsWithRes::hash
        
            typedef ::uint ( ::SireMol::ChainsWithRes::*hash_function_type)(  ) const;
            hash_function_type hash_function_value( &::SireMol::ChainsWithRes::hash );
            
            ChainsWithRes_exposer.def( 
                "hash"
                , hash_function_value
                , "Return a hash of this identifier" );
        
        }
        { //::SireMol::ChainsWithRes::isNull
        
            typedef bool ( ::SireMol::ChainsWithRes::*isNull_function_type)(  ) const;
            isNull_function_type isNull_function_value( &::SireMol::ChainsWithRes::isNull );
            
            ChainsWithRes_exposer.def( 
                "isNull"
                , isNull_function_value
                , "Is this selection null?" );
        
        }
        { //::SireMol::ChainsWithRes::map
        
            typedef ::QList< SireMol::ChainIdx > ( ::SireMol::ChainsWithRes::*map_function_type)( ::SireMol::MolInfo const & ) const;
            map_function_type map_function_value( &::SireMol::ChainsWithRes::map );
            
            ChainsWithRes_exposer.def( 
                "map"
                , map_function_value
                , ( bp::arg("molinfo") )
                , "Map this ID to the list of indicies of chains that match this ID\nThrow: SireMol::missing_chain\nThrow: SireError::invalid_index\n" );
        
        }
        ChainsWithRes_exposer.def( bp::self != bp::self );
        { //::SireMol::ChainsWithRes::operator=
        
            typedef ::SireMol::ChainsWithRes & ( ::SireMol::ChainsWithRes::*assign_function_type)( ::SireMol::ChainsWithRes const & ) ;
            assign_function_type assign_function_value( &::SireMol::ChainsWithRes::operator= );
            
            ChainsWithRes_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        ChainsWithRes_exposer.def( bp::self == bp::other< SireID::ID >() );
        ChainsWithRes_exposer.def( bp::self == bp::self );
        { //::SireMol::ChainsWithRes::resID
        
            typedef ::SireMol::ResID const & ( ::SireMol::ChainsWithRes::*resID_function_type)(  ) const;
            resID_function_type resID_function_value( &::SireMol::ChainsWithRes::resID );
            
            ChainsWithRes_exposer.def( 
                "resID"
                , resID_function_value
                , bp::return_value_policy<bp::clone_const_reference>()
                , "Return the residue ID" );
        
        }
        { //::SireMol::ChainsWithRes::toString
        
            typedef ::QString ( ::SireMol::ChainsWithRes::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::ChainsWithRes::toString );
            
            ChainsWithRes_exposer.def( 
                "toString"
                , toString_function_value
                , "Return a string representatio of this ID" );
        
        }
        { //::SireMol::ChainsWithRes::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::ChainsWithRes::typeName );
            
            ChainsWithRes_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireMol::ChainsWithRes::what
        
            typedef char const * ( ::SireMol::ChainsWithRes::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMol::ChainsWithRes::what );
            
            ChainsWithRes_exposer.def( 
                "what"
                , what_function_value
                , "" );
        
        }
        ChainsWithRes_exposer.staticmethod( "typeName" );
        ChainsWithRes_exposer.def( "__copy__", &__copy__);
        ChainsWithRes_exposer.def( "__deepcopy__", &__copy__);
        ChainsWithRes_exposer.def( "clone", &__copy__);
        ChainsWithRes_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::ChainsWithRes >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ChainsWithRes_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::ChainsWithRes >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ChainsWithRes_exposer.def( "__setstate__", &__setstate__base64< ::SireMol::ChainsWithRes > );
        ChainsWithRes_exposer.def( "__getstate__", &__getstate__base64< ::SireMol::ChainsWithRes > );
        ChainsWithRes_exposer.def( "__str__", &__str__< ::SireMol::ChainsWithRes > );
        ChainsWithRes_exposer.def( "__repr__", &__str__< ::SireMol::ChainsWithRes > );
        ChainsWithRes_exposer.def( "__hash__", &::SireMol::ChainsWithRes::hash );
    }

}
