// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "SegPropertyProperty.pypp.hpp"

namespace bp = boost::python;

#include "segproperty.hpp"

#include "segproperty.hpp"

#include "SireMaths/vector.h"

#include "SireMol/moleculeview.h"

#include "SireMol/atomidxmapping.h"

SireMol::SegProperty<SireBase::PropPtr<SireBase::Property> > __copy__(const SireMol::SegProperty<SireBase::PropPtr<SireBase::Property> > &other){ return SireMol::SegProperty<SireBase::PropPtr<SireBase::Property> >(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_SegPropertyProperty_class(){

    { //::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >
        typedef bp::class_< SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >, bp::bases< SireMol::SegProp, SireMol::MolViewProperty, SireBase::Property > > SegPropertyProperty_exposer_t;
        SegPropertyProperty_exposer_t SegPropertyProperty_exposer = SegPropertyProperty_exposer_t( "SegPropertyProperty", "", bp::init< >("") );
        bp::scope SegPropertyProperty_scope( SegPropertyProperty_exposer );
        SegPropertyProperty_exposer.def( bp::init< SireMol::MoleculeInfoData const & >(( bp::arg("molinfo") ), "") );
        SegPropertyProperty_exposer.def( bp::init< QVector< SireBase::PropPtr< SireBase::Property > > const & >(( bp::arg("values") ), "") );
        SegPropertyProperty_exposer.def( bp::init< SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > > const & >(( bp::arg("other") ), "") );
        { //::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::array
        
            typedef SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef ::QVector< SireBase::PropPtr< SireBase::Property > > const & ( ::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::*array_function_type)(  ) const;
            array_function_type array_function_value( &::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::array );
            
            SegPropertyProperty_exposer.def( 
                "array"
                , array_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::assertCanConvert
        
            typedef SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef void ( ::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::*assertCanConvert_function_type)( ::QVariant const & ) const;
            assertCanConvert_function_type assertCanConvert_function_value( &::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::assertCanConvert );
            
            SegPropertyProperty_exposer.def( 
                "assertCanConvert"
                , assertCanConvert_function_value
                , ( bp::arg("value") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::assignFrom
        
            typedef SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef void ( ::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::*assignFrom_function_type)( ::SireMol::SegProperty< QVariant > const & ) ;
            assignFrom_function_type assignFrom_function_value( &::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::assignFrom );
            
            SegPropertyProperty_exposer.def( 
                "assignFrom"
                , assignFrom_function_value
                , ( bp::arg("values") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::at
        
            typedef SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef ::SireBase::PropPtr< SireBase::Property > const & ( ::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::*at_function_type)( ::SireMol::SegIdx const & ) const;
            at_function_type at_function_value( &::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::at );
            
            SegPropertyProperty_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("segidx") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::at
        
            typedef SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef ::SireBase::PropPtr< SireBase::Property > const & ( ::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::*at_function_type)( int ) const;
            at_function_type at_function_value( &::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::at );
            
            SegPropertyProperty_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::canConvert
        
            typedef SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef bool ( ::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::*canConvert_function_type)( ::QVariant const & ) const;
            canConvert_function_type canConvert_function_value( &::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::canConvert );
            
            SegPropertyProperty_exposer.def( 
                "canConvert"
                , canConvert_function_value
                , ( bp::arg("value") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::count
        
            typedef SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef int ( ::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::*count_function_type)(  ) const;
            count_function_type count_function_value( &::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::count );
            
            SegPropertyProperty_exposer.def( 
                "count"
                , count_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::fromVariant
        
            typedef SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef ::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > > ( *fromVariant_function_type )( ::SireMol::SegProperty< QVariant > const & );
            fromVariant_function_type fromVariant_function_value( &::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::fromVariant );
            
            SegPropertyProperty_exposer.def( 
                "fromVariant"
                , fromVariant_function_value
                , ( bp::arg("values") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::get
        
            typedef SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef ::SireBase::PropPtr< SireBase::Property > const & ( ::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::*get_function_type)( ::SireMol::SegIdx const & ) const;
            get_function_type get_function_value( &::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::get );
            
            SegPropertyProperty_exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("segidx") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::get
        
            typedef SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef ::SireBase::PropPtr< SireBase::Property > const & ( ::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::*get_function_type)( int ) const;
            get_function_type get_function_value( &::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::get );
            
            SegPropertyProperty_exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::getAsProperty
        
            typedef SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef ::SireBase::PropertyPtr ( ::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::*getAsProperty_function_type)( ::SireMol::SegIdx const & ) const;
            getAsProperty_function_type getAsProperty_function_value( &::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::getAsProperty );
            
            SegPropertyProperty_exposer.def( 
                "getAsProperty"
                , getAsProperty_function_value
                , ( bp::arg("idx") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::getAsVariant
        
            typedef SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef ::QVariant ( ::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::*getAsVariant_function_type)( ::SireMol::SegIdx const & ) const;
            getAsVariant_function_type getAsVariant_function_value( &::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::getAsVariant );
            
            SegPropertyProperty_exposer.def( 
                "getAsVariant"
                , getAsVariant_function_value
                , ( bp::arg("idx") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::isCompatibleWith
        
            typedef SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef bool ( ::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::*isCompatibleWith_function_type)( ::SireMol::MoleculeInfoData const & ) const;
            isCompatibleWith_function_type isCompatibleWith_function_value( &::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::isCompatibleWith );
            
            SegPropertyProperty_exposer.def( 
                "isCompatibleWith"
                , isCompatibleWith_function_value
                , ( bp::arg("molinfo") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::isEmpty
        
            typedef SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef bool ( ::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::isEmpty );
            
            SegPropertyProperty_exposer.def( 
                "isEmpty"
                , isEmpty_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::merge
        
            typedef SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef ::SireBase::PropertyList ( ::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::*merge_function_type)( ::SireMol::MolViewProperty const &,::SireMol::AtomIdxMapping const &,::QString const &,::SireBase::PropertyMap const & ) const;
            merge_function_type merge_function_value( &::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::merge );
            
            SegPropertyProperty_exposer.def( 
                "merge"
                , merge_function_value
                , ( bp::arg("other"), bp::arg("mapping"), bp::arg("ghost")=::QString( ), bp::arg("map")=SireBase::PropertyMap() )
                , "" );
        
        }
        { //::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::nSegments
        
            typedef SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef int ( ::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::*nSegments_function_type)(  ) const;
            nSegments_function_type nSegments_function_value( &::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::nSegments );
            
            SegPropertyProperty_exposer.def( 
                "nSegments"
                , nSegments_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        SegPropertyProperty_exposer.def( bp::self != bp::self );
        { //::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::operator=
        
            typedef SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef ::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > > & ( ::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::*assign_function_type)( ::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > > const & ) ;
            assign_function_type assign_function_value( &::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::operator= );
            
            SegPropertyProperty_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        SegPropertyProperty_exposer.def( bp::self == bp::self );
        { //::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::operator[]
        
            typedef SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef ::SireBase::PropPtr< SireBase::Property > const & ( ::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::*__getitem___function_type)( ::SireMol::SegIdx const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::operator[] );
            
            SegPropertyProperty_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("segidx") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::operator[]
        
            typedef SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef ::SireBase::PropPtr< SireBase::Property > const & ( ::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::*__getitem___function_type)( int ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::operator[] );
            
            SegPropertyProperty_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("i") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::operator[]
        
            typedef SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef ::QList< SireBase::PropPtr< SireBase::Property > > ( ::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::*__getitem___function_type)( ::QList< long long > const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::operator[] );
            
            SegPropertyProperty_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("idxs") )
                , "" );
        
        }
        { //::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::operator[]
        
            typedef SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef ::QList< SireBase::PropPtr< SireBase::Property > > ( ::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::*__getitem___function_type)( ::SireBase::Slice const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::operator[] );
            
            SegPropertyProperty_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("slice") )
                , "" );
        
        }
        { //::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::set
        
            typedef SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef ::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > > & ( ::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::*set_function_type)( ::SireMol::SegIdx,::SireBase::PropPtr< SireBase::Property > const & ) ;
            set_function_type set_function_value( &::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::set );
            
            SegPropertyProperty_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("segidx"), bp::arg("value") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::size
        
            typedef SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef int ( ::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::*size_function_type)(  ) const;
            size_function_type size_function_value( &::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::size );
            
            SegPropertyProperty_exposer.def( 
                "size"
                , size_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::toString
        
            typedef SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef ::QString ( ::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::toString );
            
            SegPropertyProperty_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::toVariant
        
            typedef SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef ::SireMol::SegProperty< QVariant > ( ::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::*toVariant_function_type)(  ) const;
            toVariant_function_type toVariant_function_value( &::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::toVariant );
            
            SegPropertyProperty_exposer.def( 
                "toVariant"
                , toVariant_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::typeName
        
            typedef SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::SegProperty< SireBase::PropPtr< SireBase::Property > >::typeName );
            
            SegPropertyProperty_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        SegPropertyProperty_exposer.staticmethod( "fromVariant" );
        SegPropertyProperty_exposer.staticmethod( "typeName" );
        SegPropertyProperty_exposer.def( "__copy__", &__copy__);
        SegPropertyProperty_exposer.def( "__deepcopy__", &__copy__);
        SegPropertyProperty_exposer.def( "clone", &__copy__);
        SegPropertyProperty_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::SegProperty<SireBase::PropPtr<SireBase::Property> > >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        SegPropertyProperty_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::SegProperty<SireBase::PropPtr<SireBase::Property> > >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        SegPropertyProperty_exposer.def_pickle(sire_pickle_suite< ::SireMol::SegProperty<SireBase::PropPtr<SireBase::Property> > >());
        SegPropertyProperty_exposer.def( "__str__", &__str__< ::SireMol::SegProperty<SireBase::PropPtr<SireBase::Property> > > );
        SegPropertyProperty_exposer.def( "__repr__", &__str__< ::SireMol::SegProperty<SireBase::PropPtr<SireBase::Property> > > );
        SegPropertyProperty_exposer.def( "__len__", &__len_size< ::SireMol::SegProperty<SireBase::PropPtr<SireBase::Property> > > );
    }

}
