// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "AtomEnergies.pypp.hpp"

namespace bp = boost::python;

#include "atomenergies.h"

#include "atomenergies.h"

#include "SireMaths/vector.h"

SireMol::AtomProperty<SireUnits::Dimension::PhysUnit<1, 2, -2, 0, 0, -1, 0> > __copy__(const SireMol::AtomProperty<SireUnits::Dimension::PhysUnit<1, 2, -2, 0, 0, -1, 0> > &other){ return SireMol::AtomProperty<SireUnits::Dimension::PhysUnit<1, 2, -2, 0, 0, -1, 0> >(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/len.hpp"

void register_AtomEnergies_class(){

    { //::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >
        typedef bp::class_< SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >, bp::bases< SireMol::AtomProp, SireMol::MolViewProperty, SireBase::Property > > AtomEnergies_exposer_t;
        AtomEnergies_exposer_t AtomEnergies_exposer = AtomEnergies_exposer_t( "AtomEnergies", bp::init< >() );
        bp::scope AtomEnergies_scope( AtomEnergies_exposer );
        AtomEnergies_exposer.def( bp::init< SireMol::MoleculeInfoData const & >(( bp::arg("molinfo") )) );
        AtomEnergies_exposer.def( bp::init< SireMol::MoleculeInfoData const &, SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > const & >(( bp::arg("molinfo"), bp::arg("default_value") )) );
        AtomEnergies_exposer.def( bp::init< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > const & >(( bp::arg("value") )) );
        AtomEnergies_exposer.def( bp::init< SireBase::PackedArray2D< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > > const & >(( bp::arg("values") )) );
        AtomEnergies_exposer.def( bp::init< SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > > const & >(( bp::arg("other") )) );
        { //::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::array
        
            typedef SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > > exported_class_t;
            typedef ::SireBase::PackedArray2D< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > > const & ( ::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::*array_function_type)(  ) const;
            array_function_type array_function_value( &::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::array );
            
            AtomEnergies_exposer.def( 
                "array"
                , array_function_value
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::assertCanConvert
        
            typedef SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > > exported_class_t;
            typedef void ( ::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::*assertCanConvert_function_type)( ::QVariant const & ) const;
            assertCanConvert_function_type assertCanConvert_function_value( &::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::assertCanConvert );
            
            AtomEnergies_exposer.def( 
                "assertCanConvert"
                , assertCanConvert_function_value
                , ( bp::arg("value") ) );
        
        }
        { //::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::assignFrom
        
            typedef SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > > exported_class_t;
            typedef void ( ::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::*assignFrom_function_type)( ::SireMol::AtomProperty< QVariant > const & ) ;
            assignFrom_function_type assignFrom_function_value( &::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::assignFrom );
            
            AtomEnergies_exposer.def( 
                "assignFrom"
                , assignFrom_function_value
                , ( bp::arg("values") ) );
        
        }
        { //::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::at
        
            typedef SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > > exported_class_t;
            typedef ::SireBase::detail::PackedArray2D_Array< SireUnits::Dimension::PhysUnit< 1, 2, -0x00000000000000002, 0, 0, -0x00000000000000001, 0 > > const & ( ::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::*at_function_type)( ::SireMol::CGIdx ) const;
            at_function_type at_function_value( &::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::at );
            
            AtomEnergies_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("cgidx") )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::at
        
            typedef SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > > exported_class_t;
            typedef ::SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > const & ( ::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::*at_function_type)( ::SireMol::CGAtomIdx const & ) const;
            at_function_type at_function_value( &::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::at );
            
            AtomEnergies_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("cgatomidx") )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::canConvert
        
            typedef SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > > exported_class_t;
            typedef bool ( ::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::*canConvert_function_type)( ::QVariant const & ) const;
            canConvert_function_type canConvert_function_value( &::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::canConvert );
            
            AtomEnergies_exposer.def( 
                "canConvert"
                , canConvert_function_value
                , ( bp::arg("value") ) );
        
        }
        { //::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::copyFrom
        
            typedef SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > > exported_class_t;
            typedef void ( ::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::*copyFrom_function_type)( ::QVector< SireUnits::Dimension::PhysUnit< 1, 2, -0x00000000000000002, 0, 0, -0x00000000000000001, 0 > > const & ) ;
            copyFrom_function_type copyFrom_function_value( &::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::copyFrom );
            
            AtomEnergies_exposer.def( 
                "copyFrom"
                , copyFrom_function_value
                , ( bp::arg("values") ) );
        
        }
        { //::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::copyFrom
        
            typedef SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > > exported_class_t;
            typedef void ( ::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::*copyFrom_function_type)( ::QVector< SireUnits::Dimension::PhysUnit< 1, 2, -0x00000000000000002, 0, 0, -0x00000000000000001, 0 > > const &,::SireMol::AtomSelection const & ) ;
            copyFrom_function_type copyFrom_function_value( &::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::copyFrom );
            
            AtomEnergies_exposer.def( 
                "copyFrom"
                , copyFrom_function_value
                , ( bp::arg("values"), bp::arg("selection") ) );
        
        }
        { //::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::count
        
            typedef SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > > exported_class_t;
            typedef int ( ::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::*count_function_type)(  ) const;
            count_function_type count_function_value( &::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::count );
            
            AtomEnergies_exposer.def( 
                "count"
                , count_function_value );
        
        }
        { //::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::divide
        
            typedef SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > > exported_class_t;
            typedef ::SireBase::PropPtr< SireBase::Property > ( ::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::*divide_function_type)( ::QVector< SireMol::AtomSelection > const & ) const;
            divide_function_type divide_function_value( &::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::divide );
            
            AtomEnergies_exposer.def( 
                "divide"
                , divide_function_value
                , ( bp::arg("beads") ) );
        
        }
        { //::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::divideByResidue
        
            typedef SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > > exported_class_t;
            typedef ::SireBase::PropPtr< SireBase::Property > ( ::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::*divideByResidue_function_type)( ::SireMol::MoleculeInfoData const & ) const;
            divideByResidue_function_type divideByResidue_function_value( &::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::divideByResidue );
            
            AtomEnergies_exposer.def( 
                "divideByResidue"
                , divideByResidue_function_value
                , ( bp::arg("molinfo") ) );
        
        }
        { //::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::fromVariant
        
            typedef SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > > exported_class_t;
            typedef ::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > > ( *fromVariant_function_type )( ::SireMol::AtomProperty< QVariant > const & );
            fromVariant_function_type fromVariant_function_value( &::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::fromVariant );
            
            AtomEnergies_exposer.def( 
                "fromVariant"
                , fromVariant_function_value
                , ( bp::arg("variant") ) );
        
        }
        { //::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::get
        
            typedef SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > > exported_class_t;
            typedef ::SireBase::detail::PackedArray2D_Array< SireUnits::Dimension::PhysUnit< 1, 2, -0x00000000000000002, 0, 0, -0x00000000000000001, 0 > > const & ( ::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::*get_function_type)( ::SireMol::CGIdx ) const;
            get_function_type get_function_value( &::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::get );
            
            AtomEnergies_exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("cgidx") )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::get
        
            typedef SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > > exported_class_t;
            typedef ::SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > const & ( ::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::*get_function_type)( ::SireMol::CGAtomIdx const & ) const;
            get_function_type get_function_value( &::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::get );
            
            AtomEnergies_exposer.def( 
                "get"
                , get_function_value
                , ( bp::arg("cgatomidx") )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::isCompatibleWith
        
            typedef SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > > exported_class_t;
            typedef bool ( ::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::*isCompatibleWith_function_type)( ::SireMol::MoleculeInfoData const & ) const;
            isCompatibleWith_function_type isCompatibleWith_function_value( &::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::isCompatibleWith );
            
            AtomEnergies_exposer.def( 
                "isCompatibleWith"
                , isCompatibleWith_function_value
                , ( bp::arg("molinfo") ) );
        
        }
        { //::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::isEmpty
        
            typedef SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > > exported_class_t;
            typedef bool ( ::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::isEmpty );
            
            AtomEnergies_exposer.def( 
                "isEmpty"
                , isEmpty_function_value );
        
        }
        { //::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::matchToSelection
        
            typedef SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > > exported_class_t;
            typedef ::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > > ( ::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::*matchToSelection_function_type)( ::SireMol::AtomSelection const & ) const;
            matchToSelection_function_type matchToSelection_function_value( &::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::matchToSelection );
            
            AtomEnergies_exposer.def( 
                "matchToSelection"
                , matchToSelection_function_value
                , ( bp::arg("selection") ) );
        
        }
        { //::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::merge
        
            typedef SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > > exported_class_t;
            typedef ::SireBase::PropPtr< SireBase::Property > ( ::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::*merge_function_type)( ::SireMol::MoleculeInfoData const & ) const;
            merge_function_type merge_function_value( &::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::merge );
            
            AtomEnergies_exposer.def( 
                "merge"
                , merge_function_value
                , ( bp::arg("molinfo") ) );
        
        }
        { //::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::nAtoms
        
            typedef SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > > exported_class_t;
            typedef int ( ::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::*nAtoms_function_type)(  ) const;
            nAtoms_function_type nAtoms_function_value( &::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::nAtoms );
            
            AtomEnergies_exposer.def( 
                "nAtoms"
                , nAtoms_function_value );
        
        }
        { //::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::nAtoms
        
            typedef SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > > exported_class_t;
            typedef int ( ::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::*nAtoms_function_type)( ::SireMol::CGIdx ) const;
            nAtoms_function_type nAtoms_function_value( &::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::nAtoms );
            
            AtomEnergies_exposer.def( 
                "nAtoms"
                , nAtoms_function_value
                , ( bp::arg("cgidx") ) );
        
        }
        { //::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::nCutGroups
        
            typedef SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > > exported_class_t;
            typedef int ( ::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::*nCutGroups_function_type)(  ) const;
            nCutGroups_function_type nCutGroups_function_value( &::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::nCutGroups );
            
            AtomEnergies_exposer.def( 
                "nCutGroups"
                , nCutGroups_function_value );
        
        }
        AtomEnergies_exposer.def( bp::self != bp::self );
        { //::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::operator=
        
            typedef SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > > exported_class_t;
            typedef ::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > > & ( ::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::*assign_function_type)( ::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > > const & ) ;
            assign_function_type assign_function_value( &::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::operator= );
            
            AtomEnergies_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >() );
        
        }
        AtomEnergies_exposer.def( bp::self == bp::self );
        { //::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::operator[]
        
            typedef SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > > exported_class_t;
            typedef ::SireBase::detail::PackedArray2D_Array< SireUnits::Dimension::PhysUnit< 1, 2, -0x00000000000000002, 0, 0, -0x00000000000000001, 0 > > const & ( ::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::*__getitem___function_type)( ::SireMol::CGIdx ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::operator[] );
            
            AtomEnergies_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("cgidx") )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::operator[]
        
            typedef SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > > exported_class_t;
            typedef ::SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > const & ( ::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::*__getitem___function_type)( ::SireMol::CGAtomIdx const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::operator[] );
            
            AtomEnergies_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("cgatomidx") )
                , bp::return_value_policy< bp::copy_const_reference >() );
        
        }
        { //::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::set
        
            typedef SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > > exported_class_t;
            typedef ::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > > & ( ::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::*set_function_type)( ::SireMol::CGAtomIdx const &,::SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > const & ) ;
            set_function_type set_function_value( &::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::set );
            
            AtomEnergies_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("cgatomidx"), bp::arg("value") )
                , bp::return_self< >() );
        
        }
        { //::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::set
        
            typedef SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > > exported_class_t;
            typedef ::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > > & ( ::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::*set_function_type)( ::SireMol::CGIdx,::QVector< SireUnits::Dimension::PhysUnit< 1, 2, -0x00000000000000002, 0, 0, -0x00000000000000001, 0 > > const & ) ;
            set_function_type set_function_value( &::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::set );
            
            AtomEnergies_exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("cgidx"), bp::arg("values") )
                , bp::return_self< >() );
        
        }
        { //::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::size
        
            typedef SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > > exported_class_t;
            typedef int ( ::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::*size_function_type)(  ) const;
            size_function_type size_function_value( &::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::size );
            
            AtomEnergies_exposer.def( 
                "size"
                , size_function_value );
        
        }
        { //::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::toString
        
            typedef SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > > exported_class_t;
            typedef ::QString ( ::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::toString );
            
            AtomEnergies_exposer.def( 
                "toString"
                , toString_function_value );
        
        }
        { //::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::toVariant
        
            typedef SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > > exported_class_t;
            typedef ::SireMol::AtomProperty< QVariant > ( ::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::*toVariant_function_type)(  ) const;
            toVariant_function_type toVariant_function_value( &::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::toVariant );
            
            AtomEnergies_exposer.def( 
                "toVariant"
                , toVariant_function_value );
        
        }
        { //::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::toVector
        
            typedef SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > > exported_class_t;
            typedef ::QVector< SireUnits::Dimension::PhysUnit< 1, 2, -0x00000000000000002, 0, 0, -0x00000000000000001, 0 > > ( ::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::*toVector_function_type)(  ) const;
            toVector_function_type toVector_function_value( &::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::toVector );
            
            AtomEnergies_exposer.def( 
                "toVector"
                , toVector_function_value );
        
        }
        { //::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::toVector
        
            typedef SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > > exported_class_t;
            typedef ::QVector< SireUnits::Dimension::PhysUnit< 1, 2, -0x00000000000000002, 0, 0, -0x00000000000000001, 0 > > ( ::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::*toVector_function_type)( ::SireMol::AtomSelection const & ) const;
            toVector_function_type toVector_function_value( &::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::toVector );
            
            AtomEnergies_exposer.def( 
                "toVector"
                , toVector_function_value
                , ( bp::arg("selection") ) );
        
        }
        { //::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::typeName
        
            typedef SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::AtomProperty< SireUnits::Dimension::PhysUnit< 1, 2, -2, 0, 0, -1, 0 > >::typeName );
            
            AtomEnergies_exposer.def( 
                "typeName"
                , typeName_function_value );
        
        }
        AtomEnergies_exposer.staticmethod( "fromVariant" );
        AtomEnergies_exposer.staticmethod( "typeName" );
        AtomEnergies_exposer.def( "__copy__", &__copy__);
        AtomEnergies_exposer.def( "__deepcopy__", &__copy__);
        AtomEnergies_exposer.def( "clone", &__copy__);
        AtomEnergies_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::AtomProperty<SireUnits::Dimension::PhysUnit<1, 2, -2, 0, 0, -1, 0> > >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AtomEnergies_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::AtomProperty<SireUnits::Dimension::PhysUnit<1, 2, -2, 0, 0, -1, 0> > >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        AtomEnergies_exposer.def( "__str__", &__str__< ::SireMol::AtomProperty<SireUnits::Dimension::PhysUnit<1, 2, -2, 0, 0, -1, 0> > > );
        AtomEnergies_exposer.def( "__repr__", &__str__< ::SireMol::AtomProperty<SireUnits::Dimension::PhysUnit<1, 2, -2, 0, 0, -1, 0> > > );
        AtomEnergies_exposer.def( "__len__", &__len_size< ::SireMol::AtomProperty<SireUnits::Dimension::PhysUnit<1, 2, -2, 0, 0, -1, 0> > > );
    }

}
