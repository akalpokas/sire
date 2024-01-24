// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "Mover_Selector_Residue_.pypp.hpp"

namespace bp = boost::python;

#include "mover_metaid.h"

#include "mover_metaid.h"

#include "SireMaths/align.h"

#include "SireMaths/axisset.h"

#include "SireMaths/matrix.h"

#include "SireMaths/quaternion.h"

#include "SireMaths/rotate.h"

#include "SireMaths/vectorproperty.h"

#include "SireMol/errors.h"

#include "SireUnits/units.h"

#include "SireVol/coordgroup.h"

#include "SireVol/space.h"

#include "angleid.h"

#include "atomcoords.h"

#include "atommatcher.h"

#include "atommatchers.h"

#include "atomvelocities.h"

#include "bondid.h"

#include "connectivity.h"

#include "dihedralid.h"

#include "improperid.h"

#include "mover.h"

#include "tostring.h"

#include "weightfunction.h"

#include "mover.h"

SireMol::Mover<SireMol::Selector<SireMol::Residue> > __copy__(const SireMol::Mover<SireMol::Selector<SireMol::Residue> > &other){ return SireMol::Mover<SireMol::Selector<SireMol::Residue> >(other); }

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_Mover_Selector_Residue__class(){

    { //::SireMol::Mover< SireMol::Selector< SireMol::Residue > >
        typedef bp::class_< SireMol::Mover< SireMol::Selector< SireMol::Residue > >, bp::bases< SireMol::MoverBase, SireMol::Selector<SireMol::Residue>, SireMol::MoleculeView, SireBase::Property > > Mover_Selector_Residue__exposer_t;
        Mover_Selector_Residue__exposer_t Mover_Selector_Residue__exposer = Mover_Selector_Residue__exposer_t( "Mover_Selector_Residue_", "", bp::init< >("") );
        bp::scope Mover_Selector_Residue__scope( Mover_Selector_Residue__exposer );
        Mover_Selector_Residue__exposer.def( bp::init< SireMol::Selector< SireMol::Residue > const & >(( bp::arg("view") ), "") );
        Mover_Selector_Residue__exposer.def( bp::init< SireMol::Selector< SireMol::Residue > const &, SireMol::AtomSelection const & >(( bp::arg("view"), bp::arg("movable_atoms") ), "") );
        Mover_Selector_Residue__exposer.def( bp::init< SireMol::Mover< SireMol::Selector< SireMol::Residue > > const & >(( bp::arg("other") ), "") );
        { //::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::align
        
            typedef SireMol::Mover< SireMol::Selector< SireMol::Residue > > exported_class_t;
            typedef ::SireMol::Mover< SireMol::Selector< SireMol::Residue > > & ( ::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::*align_function_type)( ::SireMol::MoleculeView const &,::SireBase::PropertyMap const & ) ;
            align_function_type align_function_value( &::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::align );
            
            Mover_Selector_Residue__exposer.def( 
                "align"
                , align_function_value
                , ( bp::arg("other"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::align
        
            typedef SireMol::Mover< SireMol::Selector< SireMol::Residue > > exported_class_t;
            typedef ::SireMol::Mover< SireMol::Selector< SireMol::Residue > > & ( ::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::*align_function_type)( ::SireMol::MoleculeView const &,::SireBase::PropertyMap const &,::SireBase::PropertyMap const & ) ;
            align_function_type align_function_value( &::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::align );
            
            Mover_Selector_Residue__exposer.def( 
                "align"
                , align_function_value
                , ( bp::arg("other"), bp::arg("map0"), bp::arg("map1") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::align
        
            typedef SireMol::Mover< SireMol::Selector< SireMol::Residue > > exported_class_t;
            typedef ::SireMol::Mover< SireMol::Selector< SireMol::Residue > > & ( ::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::*align_function_type)( ::SireMol::MoleculeView const &,::SireMol::AtomMatcher const &,::SireBase::PropertyMap const & ) ;
            align_function_type align_function_value( &::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::align );
            
            Mover_Selector_Residue__exposer.def( 
                "align"
                , align_function_value
                , ( bp::arg("other"), bp::arg("matcher"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::align
        
            typedef SireMol::Mover< SireMol::Selector< SireMol::Residue > > exported_class_t;
            typedef ::SireMol::Mover< SireMol::Selector< SireMol::Residue > > & ( ::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::*align_function_type)( ::SireMol::MoleculeView const &,::SireMol::AtomMatcher const &,::SireBase::PropertyMap const &,::SireBase::PropertyMap const & ) ;
            align_function_type align_function_value( &::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::align );
            
            Mover_Selector_Residue__exposer.def( 
                "align"
                , align_function_value
                , ( bp::arg("other"), bp::arg("matcher"), bp::arg("map0"), bp::arg("map1") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::change
        
            typedef SireMol::Mover< SireMol::Selector< SireMol::Residue > > exported_class_t;
            typedef ::SireMol::Mover< SireMol::Selector< SireMol::Residue > > & ( ::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::*change_function_type)( ::SireMol::BondID const &,::SireUnits::Dimension::Length,::SireBase::PropertyMap const & ) ;
            change_function_type change_function_value( &::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::change );
            
            Mover_Selector_Residue__exposer.def( 
                "change"
                , change_function_value
                , ( bp::arg("bond"), bp::arg("delta"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::change
        
            typedef SireMol::Mover< SireMol::Selector< SireMol::Residue > > exported_class_t;
            typedef ::SireMol::Mover< SireMol::Selector< SireMol::Residue > > & ( ::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::*change_function_type)( ::SireMol::AngleID const &,::SireUnits::Dimension::Angle,::SireBase::PropertyMap const & ) ;
            change_function_type change_function_value( &::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::change );
            
            Mover_Selector_Residue__exposer.def( 
                "change"
                , change_function_value
                , ( bp::arg("angle"), bp::arg("delta"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::change
        
            typedef SireMol::Mover< SireMol::Selector< SireMol::Residue > > exported_class_t;
            typedef ::SireMol::Mover< SireMol::Selector< SireMol::Residue > > & ( ::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::*change_function_type)( ::SireMol::DihedralID const &,::SireUnits::Dimension::Angle,::SireBase::PropertyMap const & ) ;
            change_function_type change_function_value( &::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::change );
            
            Mover_Selector_Residue__exposer.def( 
                "change"
                , change_function_value
                , ( bp::arg("dihedral"), bp::arg("delta"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::change
        
            typedef SireMol::Mover< SireMol::Selector< SireMol::Residue > > exported_class_t;
            typedef ::SireMol::Mover< SireMol::Selector< SireMol::Residue > > & ( ::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::*change_function_type)( ::SireMol::BondID const &,::SireUnits::Dimension::Angle,::SireBase::PropertyMap const & ) ;
            change_function_type change_function_value( &::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::change );
            
            Mover_Selector_Residue__exposer.def( 
                "change"
                , change_function_value
                , ( bp::arg("bond"), bp::arg("delta"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::change
        
            typedef SireMol::Mover< SireMol::Selector< SireMol::Residue > > exported_class_t;
            typedef ::SireMol::Mover< SireMol::Selector< SireMol::Residue > > & ( ::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::*change_function_type)( ::SireMol::ImproperID const &,::SireUnits::Dimension::Angle,::SireBase::PropertyMap const & ) ;
            change_function_type change_function_value( &::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::change );
            
            Mover_Selector_Residue__exposer.def( 
                "change"
                , change_function_value
                , ( bp::arg("improper"), bp::arg("delta"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::changeFrame
        
            typedef SireMol::Mover< SireMol::Selector< SireMol::Residue > > exported_class_t;
            typedef ::SireMol::Mover< SireMol::Selector< SireMol::Residue > > & ( ::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::*changeFrame_function_type)( ::SireMaths::AxisSet const &,::SireMaths::AxisSet const &,::SireBase::PropertyMap const & ) ;
            changeFrame_function_type changeFrame_function_value( &::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::changeFrame );
            
            Mover_Selector_Residue__exposer.def( 
                "changeFrame"
                , changeFrame_function_value
                , ( bp::arg("from_frame"), bp::arg("to_frame"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::commit
        
            typedef SireMol::Mover< SireMol::Selector< SireMol::Residue > > exported_class_t;
            typedef ::SireMol::Selector< SireMol::Residue > ( ::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::*commit_function_type)(  ) const;
            commit_function_type commit_function_value( &::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::commit );
            
            Mover_Selector_Residue__exposer.def( 
                "commit"
                , commit_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::makeWhole
        
            typedef SireMol::Mover< SireMol::Selector< SireMol::Residue > > exported_class_t;
            typedef ::SireMol::Mover< SireMol::Selector< SireMol::Residue > > & ( ::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::*makeWhole_function_type)( ::SireBase::PropertyMap const & ) ;
            makeWhole_function_type makeWhole_function_value( &::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::makeWhole );
            
            Mover_Selector_Residue__exposer.def( 
                "makeWhole"
                , makeWhole_function_value
                , ( bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::makeWhole
        
            typedef SireMol::Mover< SireMol::Selector< SireMol::Residue > > exported_class_t;
            typedef ::SireMol::Mover< SireMol::Selector< SireMol::Residue > > & ( ::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::*makeWhole_function_type)( ::SireMaths::Vector const &,::SireBase::PropertyMap const & ) ;
            makeWhole_function_type makeWhole_function_value( &::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::makeWhole );
            
            Mover_Selector_Residue__exposer.def( 
                "makeWhole"
                , makeWhole_function_value
                , ( bp::arg("center"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::mapInto
        
            typedef SireMol::Mover< SireMol::Selector< SireMol::Residue > > exported_class_t;
            typedef ::SireMol::Mover< SireMol::Selector< SireMol::Residue > > & ( ::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::*mapInto_function_type)( ::SireMaths::AxisSet const &,::SireBase::PropertyMap const & ) ;
            mapInto_function_type mapInto_function_value( &::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::mapInto );
            
            Mover_Selector_Residue__exposer.def( 
                "mapInto"
                , mapInto_function_value
                , ( bp::arg("axes"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::operator=
        
            typedef SireMol::Mover< SireMol::Selector< SireMol::Residue > > exported_class_t;
            typedef ::SireMol::Mover< SireMol::Selector< SireMol::Residue > > & ( ::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::*assign_function_type)( ::SireMol::Mover< SireMol::Selector< SireMol::Residue > > const & ) ;
            assign_function_type assign_function_value( &::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::operator= );
            
            Mover_Selector_Residue__exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::operator=
        
            typedef SireMol::Mover< SireMol::Selector< SireMol::Residue > > exported_class_t;
            typedef ::SireMol::Mover< SireMol::Selector< SireMol::Residue > > & ( ::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::*assign_function_type)( ::SireMol::Selector< SireMol::Residue > const & ) ;
            assign_function_type assign_function_value( &::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::operator= );
            
            Mover_Selector_Residue__exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::rotate
        
            typedef SireMol::Mover< SireMol::Selector< SireMol::Residue > > exported_class_t;
            typedef ::SireMol::Mover< SireMol::Selector< SireMol::Residue > > & ( ::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::*rotate_function_type)( ::SireMaths::Quaternion const &,::SireMaths::Vector const &,::SireBase::PropertyMap const & ) ;
            rotate_function_type rotate_function_value( &::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::rotate );
            
            Mover_Selector_Residue__exposer.def( 
                "rotate"
                , rotate_function_value
                , ( bp::arg("quat"), bp::arg("point"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::rotate
        
            typedef SireMol::Mover< SireMol::Selector< SireMol::Residue > > exported_class_t;
            typedef ::SireMol::Mover< SireMol::Selector< SireMol::Residue > > & ( ::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::*rotate_function_type)( ::SireMaths::Matrix const &,::SireMaths::Vector const &,::SireBase::PropertyMap const & ) ;
            rotate_function_type rotate_function_value( &::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::rotate );
            
            Mover_Selector_Residue__exposer.def( 
                "rotate"
                , rotate_function_value
                , ( bp::arg("rotmat"), bp::arg("point"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::set
        
            typedef SireMol::Mover< SireMol::Selector< SireMol::Residue > > exported_class_t;
            typedef ::SireMol::Mover< SireMol::Selector< SireMol::Residue > > & ( ::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::*set_function_type)( ::SireMol::BondID const &,::SireUnits::Dimension::Length,::SireBase::PropertyMap const & ) ;
            set_function_type set_function_value( &::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::set );
            
            Mover_Selector_Residue__exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("bond"), bp::arg("value"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::set
        
            typedef SireMol::Mover< SireMol::Selector< SireMol::Residue > > exported_class_t;
            typedef ::SireMol::Mover< SireMol::Selector< SireMol::Residue > > & ( ::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::*set_function_type)( ::SireMol::AngleID const &,::SireUnits::Dimension::Angle,::SireBase::PropertyMap const & ) ;
            set_function_type set_function_value( &::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::set );
            
            Mover_Selector_Residue__exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("angle"), bp::arg("value"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::set
        
            typedef SireMol::Mover< SireMol::Selector< SireMol::Residue > > exported_class_t;
            typedef ::SireMol::Mover< SireMol::Selector< SireMol::Residue > > & ( ::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::*set_function_type)( ::SireMol::DihedralID const &,::SireUnits::Dimension::Angle,::SireBase::PropertyMap const & ) ;
            set_function_type set_function_value( &::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::set );
            
            Mover_Selector_Residue__exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("dihedral"), bp::arg("value"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::set
        
            typedef SireMol::Mover< SireMol::Selector< SireMol::Residue > > exported_class_t;
            typedef ::SireMol::Mover< SireMol::Selector< SireMol::Residue > > & ( ::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::*set_function_type)( ::SireMol::ImproperID const &,::SireUnits::Dimension::Angle,::SireBase::PropertyMap const & ) ;
            set_function_type set_function_value( &::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::set );
            
            Mover_Selector_Residue__exposer.def( 
                "set"
                , set_function_value
                , ( bp::arg("improper"), bp::arg("value"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::setAll
        
            typedef SireMol::Mover< SireMol::Selector< SireMol::Residue > > exported_class_t;
            typedef ::SireMol::Mover< SireMol::Selector< SireMol::Residue > > & ( ::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::*setAll_function_type)( ::SireMol::DihedralID const &,::SireUnits::Dimension::Angle,::SireBase::PropertyMap const & ) ;
            setAll_function_type setAll_function_value( &::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::setAll );
            
            Mover_Selector_Residue__exposer.def( 
                "setAll"
                , setAll_function_value
                , ( bp::arg("dihedral"), bp::arg("value"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::toString
        
            typedef SireMol::Mover< SireMol::Selector< SireMol::Residue > > exported_class_t;
            typedef ::QString ( ::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::toString );
            
            Mover_Selector_Residue__exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::transform
        
            typedef SireMol::Mover< SireMol::Selector< SireMol::Residue > > exported_class_t;
            typedef ::SireMol::Mover< SireMol::Selector< SireMol::Residue > > & ( ::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::*transform_function_type)( ::SireMaths::Transform const &,::SireBase::PropertyMap const & ) ;
            transform_function_type transform_function_value( &::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::transform );
            
            Mover_Selector_Residue__exposer.def( 
                "transform"
                , transform_function_value
                , ( bp::arg("transform"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::translate
        
            typedef SireMol::Mover< SireMol::Selector< SireMol::Residue > > exported_class_t;
            typedef ::SireMol::Mover< SireMol::Selector< SireMol::Residue > > & ( ::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::*translate_function_type)( ::SireMaths::Vector const &,::SireBase::PropertyMap const & ) ;
            translate_function_type translate_function_value( &::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::translate );
            
            Mover_Selector_Residue__exposer.def( 
                "translate"
                , translate_function_value
                , ( bp::arg("delta"), bp::arg("map")=SireBase::PropertyMap() )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::typeName
        
            typedef SireMol::Mover< SireMol::Selector< SireMol::Residue > > exported_class_t;
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::Mover< SireMol::Selector< SireMol::Residue > >::typeName );
            
            Mover_Selector_Residue__exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        Mover_Selector_Residue__exposer.staticmethod( "typeName" );
        Mover_Selector_Residue__exposer.def( "__copy__", &__copy__);
        Mover_Selector_Residue__exposer.def( "__deepcopy__", &__copy__);
        Mover_Selector_Residue__exposer.def( "clone", &__copy__);
        Mover_Selector_Residue__exposer.def( "__str__", &__str__< ::SireMol::Mover<SireMol::Selector<SireMol::Residue> > > );
        Mover_Selector_Residue__exposer.def( "__repr__", &__str__< ::SireMol::Mover<SireMol::Selector<SireMol::Residue> > > );
        Mover_Selector_Residue__exposer.def( "__len__", &__len_size< ::SireMol::Mover<SireMol::Selector<SireMol::Residue> > > );
    }

}
