// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "ResStringProperty.pypp.hpp"

namespace bp = boost::python;

#include "resproperty.hpp"

#include "resproperty.hpp"

#include "SireMaths/vector.h"

#include "SireMol/moleculeview.h"

SireMol::ResProperty<QString> __copy__(const SireMol::ResProperty<QString> &other){ return SireMol::ResProperty<QString>(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_ResStringProperty_class(){

    { //::SireMol::ResProperty< QString >
        typedef bp::class_< SireMol::ResProperty< QString >, bp::bases< SireMol::ResProp, SireMol::MolViewProperty, SireBase::Property > > ResStringProperty_exposer_t;
        ResStringProperty_exposer_t ResStringProperty_exposer = ResStringProperty_exposer_t( "ResStringProperty", "", bp::init< >("") );
        bp::scope ResStringProperty_scope( ResStringProperty_exposer );
        ResStringProperty_exposer.def( bp::init< SireMol::MoleculeInfoData const & >(( bp::arg("molinfo") ), "") );
        ResStringProperty_exposer.def( bp::init< QVector< QString > const & >(( bp::arg("values") ), "") );
        ResStringProperty_exposer.def( bp::init< SireMol::ResProperty< QString > const & >(( bp::arg("other") ), "") );
        { //::SireMol::ResProperty< QString >::array
        
            typedef SireMol::ResProperty< QString > exported_class_t;
            typedef ::QVector< QString > const & ( ::SireMol::ResProperty< QString >::*array_function_type)(  ) const;
            array_function_type array_function_value( &::SireMol::ResProperty< QString >::array );
            
            ResStringProperty_exposer.def( 
                "array"
                , array_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::ResProperty< QString >::assertCanConvert
        
            typedef SireMol::ResProperty< QString > exported_class_t;
            typedef void ( ::SireMol::ResProperty< QString >::*assertCanConvert_function_type)( ::QVariant const & ) const;
            assertCanConvert_function_type assertCanConvert_function_value( &::SireMol::ResProperty< QString >::assertCanConvert );
            
            ResStringProperty_exposer.def( 
                "assertCanConvert"
                , assertCanConvert_function_value
                , ( bp::arg("value") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ResProperty< QString >::assignFrom
        
            typedef SireMol::ResProperty< QString > exported_class_t;
            typedef void ( ::SireMol::ResProperty< QString >::*assignFrom_function_type)( ::SireMol::ResProperty< QVariant > const & ) ;
            assignFrom_function_type assignFrom_function_value( &::SireMol::ResProperty< QString >::assignFrom );
            
            ResStringProperty_exposer.def( 
                "assignFrom"
                , assignFrom_function_value
                , ( bp::arg("values") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ResProperty< QString >::at
        
            typedef SireMol::ResProperty< QString > exported_class_t;
            typedef ::QString const & ( ::SireMol::ResProperty< QString >::*at_function_type)( ::SireMol::ResIdx const & ) const;
            at_function_type at_function_value( &::SireMol::ResProperty< QString >::at );
            
            ResStringProperty_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("residx") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::ResProperty< QString >::at
        
            typedef SireMol::ResProperty< QString > exported_class_t;
            typedef ::QString const & ( ::SireMol::ResProperty< QString >::*at_function_type)( int ) const;
            at_function_type at_function_value( &::SireMol::ResProperty< QString >::at );
            
            ResStringProperty_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::ResProperty< QString >::canConvert
        
            typedef SireMol::ResProperty< QString > exported_class_t;
            typedef bool ( ::SireMol::ResProperty< QString >::*canConvert_function_type)( ::QVariant const & ) const;
            canConvert_function_type canConvert_function_value( &::SireMol::ResProperty< QString >::canConvert );
            
            ResStringProperty_exposer.def( 
                "canConvert"
                , canConvert_function_value
                , ( bp::arg("value") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ResProperty< QString >::count
        
            typedef SireMol::ResProperty< QString > exported_class_t;
            typedef int ( ::SireMol::ResProperty< QString >::*count_function_type)(  ) const;
            count_function_type count_function_value( &::SireMol::ResProperty< QString >::count );
            
            ResStringProperty_exposer.def( 
                "count"
                , count_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ResProperty< QString >::fromVariant
        
            typedef SireMol::ResProperty< QString > exported_class_t;
            typedef ::SireMol::ResProperty< QString > ( *fromVariant_function_type )( ::SireMol::ResProperty< QVariant > const & );
            fromVariant_function_type fromVariant_function_value( &::SireMol::ResProperty< QString >::fromVariant );
            
            ResStringProperty_exposer.def( 
                "fromVariant"
                , fromVariant_function_value
                , ( bp::arg("values") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ResProperty< QString >::get
        
            typedef SireMol::ResProperty< QString > exported_class_t;
            typedef ::QString const & ( ::SireMol::ResProperty< QString >::*get_function_type)( ::SireMol::ResIdx const & ) const;
            get_function_type get_function_value( &::SireMol::ResProperty< QString >::get );
            
            ResStringProperty_exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("residx") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::ResProperty< QString >::get
        
            typedef SireMol::ResProperty< QString > exported_class_t;
            typedef ::QString const & ( ::SireMol::ResProperty< QString >::*get_function_type)( int ) const;
            get_function_type get_function_value( &::SireMol::ResProperty< QString >::get );
            
            ResStringProperty_exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::ResProperty< QString >::getAsProperty
        
            typedef SireMol::ResProperty< QString > exported_class_t;
            typedef ::SireBase::PropertyPtr ( ::SireMol::ResProperty< QString >::*getAsProperty_function_type)( ::SireMol::ResIdx const & ) const;
            getAsProperty_function_type getAsProperty_function_value( &::SireMol::ResProperty< QString >::getAsProperty );
            
            ResStringProperty_exposer.def( 
                "getAsProperty"
                , getAsProperty_function_value
                , ( bp::arg("idx") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ResProperty< QString >::getAsVariant
        
            typedef SireMol::ResProperty< QString > exported_class_t;
            typedef ::QVariant ( ::SireMol::ResProperty< QString >::*getAsVariant_function_type)( ::SireMol::ResIdx const & ) const;
            getAsVariant_function_type getAsVariant_function_value( &::SireMol::ResProperty< QString >::getAsVariant );
            
            ResStringProperty_exposer.def( 
                "getAsVariant"
                , getAsVariant_function_value
                , ( bp::arg("idx") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ResProperty< QString >::isCompatibleWith
        
            typedef SireMol::ResProperty< QString > exported_class_t;
            typedef bool ( ::SireMol::ResProperty< QString >::*isCompatibleWith_function_type)( ::SireMol::MoleculeInfoData const & ) const;
            isCompatibleWith_function_type isCompatibleWith_function_value( &::SireMol::ResProperty< QString >::isCompatibleWith );
            
            ResStringProperty_exposer.def( 
                "isCompatibleWith"
                , isCompatibleWith_function_value
                , ( bp::arg("molinfo") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ResProperty< QString >::isEmpty
        
            typedef SireMol::ResProperty< QString > exported_class_t;
            typedef bool ( ::SireMol::ResProperty< QString >::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMol::ResProperty< QString >::isEmpty );
            
            ResStringProperty_exposer.def( 
                "isEmpty"
                , isEmpty_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ResProperty< QString >::nResidues
        
            typedef SireMol::ResProperty< QString > exported_class_t;
            typedef int ( ::SireMol::ResProperty< QString >::*nResidues_function_type)(  ) const;
            nResidues_function_type nResidues_function_value( &::SireMol::ResProperty< QString >::nResidues );
            
            ResStringProperty_exposer.def( 
                "nResidues"
                , nResidues_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        ResStringProperty_exposer.def( bp::self != bp::self );
        { //::SireMol::ResProperty< QString >::operator=
        
            typedef SireMol::ResProperty< QString > exported_class_t;
            typedef ::SireMol::ResProperty< QString > & ( ::SireMol::ResProperty< QString >::*assign_function_type)( ::SireMol::ResProperty< QString > const & ) ;
            assign_function_type assign_function_value( &::SireMol::ResProperty< QString >::operator= );
            
            ResStringProperty_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        ResStringProperty_exposer.def( bp::self == bp::self );
        { //::SireMol::ResProperty< QString >::operator[]
        
            typedef SireMol::ResProperty< QString > exported_class_t;
            typedef ::QString const & ( ::SireMol::ResProperty< QString >::*__getitem___function_type)( ::SireMol::ResIdx const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::ResProperty< QString >::operator[] );
            
            ResStringProperty_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("residx") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::ResProperty< QString >::operator[]
        
            typedef SireMol::ResProperty< QString > exported_class_t;
            typedef ::QString const & ( ::SireMol::ResProperty< QString >::*__getitem___function_type)( int ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::ResProperty< QString >::operator[] );
            
            ResStringProperty_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::ResProperty< QString >::operator[]
        
            typedef SireMol::ResProperty< QString > exported_class_t;
            typedef ::QList< QString > ( ::SireMol::ResProperty< QString >::*__getitem___function_type)( ::QList< long long > const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::ResProperty< QString >::operator[] );
            
            ResStringProperty_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("idxs") )
                , "" );
        
        }
        { //::SireMol::ResProperty< QString >::operator[]
        
            typedef SireMol::ResProperty< QString > exported_class_t;
            typedef ::QList< QString > ( ::SireMol::ResProperty< QString >::*__getitem___function_type)( ::SireBase::Slice const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::ResProperty< QString >::operator[] );
            
            ResStringProperty_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("slice") )
                , "" );
        
        }
        { //::SireMol::ResProperty< QString >::set
        
            typedef SireMol::ResProperty< QString > exported_class_t;
            typedef ::SireMol::ResProperty< QString > & ( ::SireMol::ResProperty< QString >::*set_function_type)( ::SireMol::ResIdx,::QString const & ) ;
            set_function_type set_function_value( &::SireMol::ResProperty< QString >::set );
            
            ResStringProperty_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("residx"), bp::arg("value") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::ResProperty< QString >::size
        
            typedef SireMol::ResProperty< QString > exported_class_t;
            typedef int ( ::SireMol::ResProperty< QString >::*size_function_type)(  ) const;
            size_function_type size_function_value( &::SireMol::ResProperty< QString >::size );
            
            ResStringProperty_exposer.def( 
                "size"
                , size_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ResProperty< QString >::toString
        
            typedef SireMol::ResProperty< QString > exported_class_t;
            typedef ::QString ( ::SireMol::ResProperty< QString >::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::ResProperty< QString >::toString );
            
            ResStringProperty_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ResProperty< QString >::toVariant
        
            typedef SireMol::ResProperty< QString > exported_class_t;
            typedef ::SireMol::ResProperty< QVariant > ( ::SireMol::ResProperty< QString >::*toVariant_function_type)(  ) const;
            toVariant_function_type toVariant_function_value( &::SireMol::ResProperty< QString >::toVariant );
            
            ResStringProperty_exposer.def( 
                "toVariant"
                , toVariant_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::ResProperty< QString >::typeName
        
            typedef SireMol::ResProperty< QString > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::ResProperty< QString >::typeName );
            
            ResStringProperty_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        ResStringProperty_exposer.staticmethod( "fromVariant" );
        ResStringProperty_exposer.staticmethod( "typeName" );
        ResStringProperty_exposer.def( "__copy__", &__copy__);
        ResStringProperty_exposer.def( "__deepcopy__", &__copy__);
        ResStringProperty_exposer.def( "clone", &__copy__);
        ResStringProperty_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::ResProperty<QString> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ResStringProperty_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::ResProperty<QString> >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        ResStringProperty_exposer.def_pickle(sire_pickle_suite< ::SireMol::ResProperty<QString> >());
        ResStringProperty_exposer.def( "__str__", &__str__< ::SireMol::ResProperty<QString> > );
        ResStringProperty_exposer.def( "__repr__", &__str__< ::SireMol::ResProperty<QString> > );
        ResStringProperty_exposer.def( "__len__", &__len_size< ::SireMol::ResProperty<QString> > );
    }

}
