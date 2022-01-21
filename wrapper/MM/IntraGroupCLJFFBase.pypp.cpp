// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 2 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "IntraGroupCLJFFBase.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/mover.hpp"

#include "SireMol/partialmolecule.h"

#include "intracljff.h"

#include "intracljff.h"

SireFF::Intra2B2GFF<SireMM::CLJPotentialInterface<SireMM::IntraCLJPotential> > __copy__(const SireFF::Intra2B2GFF<SireMM::CLJPotentialInterface<SireMM::IntraCLJPotential> > &other){ return SireFF::Intra2B2GFF<SireMM::CLJPotentialInterface<SireMM::IntraCLJPotential> >(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/len.hpp"

void register_IntraGroupCLJFFBase_class(){

    { //::SireFF::Intra2B2GFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >
        typedef bp::class_< SireFF::Intra2B2GFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >, bp::bases< SireMM::CLJPotentialInterface<SireMM::IntraCLJPotential>, SireFF::G2FF, SireFF::FF, SireMol::MolGroupsBase, SireBase::Property > > IntraGroupCLJFFBase_exposer_t;
        IntraGroupCLJFFBase_exposer_t IntraGroupCLJFFBase_exposer = IntraGroupCLJFFBase_exposer_t( "IntraGroupCLJFFBase", "", bp::init< >("") );
        bp::scope IntraGroupCLJFFBase_scope( IntraGroupCLJFFBase_exposer );
        IntraGroupCLJFFBase_exposer.def( bp::init< QString const & >(( bp::arg("name") ), "") );
        IntraGroupCLJFFBase_exposer.def( bp::init< SireFF::Intra2B2GFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > > const & >(( bp::arg("other") ), "") );
        { //::SireFF::Intra2B2GFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::components
        
            typedef SireFF::Intra2B2GFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > > exported_class_t;
            typedef ::SireFF::Intra2B2GFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::Components const & ( ::SireFF::Intra2B2GFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::*components_function_type)(  ) const;
            components_function_type components_function_value( &::SireFF::Intra2B2GFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::components );
            
            IntraGroupCLJFFBase_exposer.def( 
                "components"
                , components_function_value
                , bp::return_value_policy<bp::clone_const_reference>()
                , "" );
        
        }
        { //::SireFF::Intra2B2GFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::containsProperty
        
            typedef SireFF::Intra2B2GFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > > exported_class_t;
            typedef bool ( ::SireFF::Intra2B2GFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::*containsProperty_function_type)( ::QString const & ) const;
            containsProperty_function_type containsProperty_function_value( &::SireFF::Intra2B2GFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::containsProperty );
            
            IntraGroupCLJFFBase_exposer.def( 
                "containsProperty"
                , containsProperty_function_value
                , ( bp::arg("name") )
                , "" );
        
        }
        { //::SireFF::Intra2B2GFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::mustNowRecalculateFromScratch
        
            typedef SireFF::Intra2B2GFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > > exported_class_t;
            typedef void ( ::SireFF::Intra2B2GFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::*mustNowRecalculateFromScratch_function_type)(  ) ;
            mustNowRecalculateFromScratch_function_type mustNowRecalculateFromScratch_function_value( &::SireFF::Intra2B2GFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::mustNowRecalculateFromScratch );
            
            IntraGroupCLJFFBase_exposer.def( 
                "mustNowRecalculateFromScratch"
                , mustNowRecalculateFromScratch_function_value
                , "" );
        
        }
        IntraGroupCLJFFBase_exposer.def( bp::self != bp::self );
        { //::SireFF::Intra2B2GFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::operator=
        
            typedef SireFF::Intra2B2GFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > > exported_class_t;
            typedef ::SireFF::Intra2B2GFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > > & ( ::SireFF::Intra2B2GFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::*assign_function_type)( ::SireFF::Intra2B2GFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > > const & ) ;
            assign_function_type assign_function_value( &::SireFF::Intra2B2GFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::operator= );
            
            IntraGroupCLJFFBase_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        IntraGroupCLJFFBase_exposer.def( bp::self == bp::self );
        { //::SireFF::Intra2B2GFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::properties
        
            typedef SireFF::Intra2B2GFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > > exported_class_t;
            typedef ::SireBase::Properties const & ( ::SireFF::Intra2B2GFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::*properties_function_type)(  ) const;
            properties_function_type properties_function_value( &::SireFF::Intra2B2GFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::properties );
            
            IntraGroupCLJFFBase_exposer.def( 
                "properties"
                , properties_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireFF::Intra2B2GFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::property
        
            typedef SireFF::Intra2B2GFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > > exported_class_t;
            typedef ::SireBase::Property const & ( ::SireFF::Intra2B2GFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::*property_function_type)( ::QString const & ) const;
            property_function_type property_function_value( &::SireFF::Intra2B2GFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::property );
            
            IntraGroupCLJFFBase_exposer.def( 
                "property"
                , property_function_value
                , ( bp::arg("name") )
                , bp::return_value_policy<bp::clone_const_reference>()
                , "" );
        
        }
        { //::SireFF::Intra2B2GFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::setProperty
        
            typedef SireFF::Intra2B2GFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > > exported_class_t;
            typedef bool ( ::SireFF::Intra2B2GFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::*setProperty_function_type)( ::QString const &,::SireBase::Property const & ) ;
            setProperty_function_type setProperty_function_value( &::SireFF::Intra2B2GFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::setProperty );
            
            IntraGroupCLJFFBase_exposer.def( 
                "setProperty"
                , setProperty_function_value
                , ( bp::arg("name"), bp::arg("property") )
                , "" );
        
        }
        { //::SireFF::Intra2B2GFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::typeName
        
            typedef SireFF::Intra2B2GFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireFF::Intra2B2GFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::typeName );
            
            IntraGroupCLJFFBase_exposer.def( 
                "typeName"
                , typeName_function_value
                , "" );
        
        }
        { //::SireFF::Intra2B2GFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::what
        
            typedef SireFF::Intra2B2GFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > > exported_class_t;
            typedef char const * ( ::SireFF::Intra2B2GFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireFF::Intra2B2GFF< SireMM::CLJPotentialInterface< SireMM::IntraCLJPotential > >::what );
            
            IntraGroupCLJFFBase_exposer.def( 
                "what"
                , what_function_value
                , "" );
        
        }
        IntraGroupCLJFFBase_exposer.staticmethod( "typeName" );
        IntraGroupCLJFFBase_exposer.def( "__copy__", &__copy__);
        IntraGroupCLJFFBase_exposer.def( "__deepcopy__", &__copy__);
        IntraGroupCLJFFBase_exposer.def( "clone", &__copy__);
        IntraGroupCLJFFBase_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireFF::Intra2B2GFF<SireMM::CLJPotentialInterface<SireMM::IntraCLJPotential> > >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        IntraGroupCLJFFBase_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireFF::Intra2B2GFF<SireMM::CLJPotentialInterface<SireMM::IntraCLJPotential> > >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        IntraGroupCLJFFBase_exposer.def( "__setstate__", &__setstate__base64< ::SireFF::Intra2B2GFF<SireMM::CLJPotentialInterface<SireMM::IntraCLJPotential> > > );
        IntraGroupCLJFFBase_exposer.def( "__getstate__", &__getstate__base64< ::SireFF::Intra2B2GFF<SireMM::CLJPotentialInterface<SireMM::IntraCLJPotential> > > );
        IntraGroupCLJFFBase_exposer.def( "__str__", &__str__< ::SireFF::Intra2B2GFF<SireMM::CLJPotentialInterface<SireMM::IntraCLJPotential> > > );
        IntraGroupCLJFFBase_exposer.def( "__repr__", &__str__< ::SireFF::Intra2B2GFF<SireMM::CLJPotentialInterface<SireMM::IntraCLJPotential> > > );
        IntraGroupCLJFFBase_exposer.def( "__len__", &__len_count< ::SireFF::Intra2B2GFF<SireMM::CLJPotentialInterface<SireMM::IntraCLJPotential> > > );
    }

}
