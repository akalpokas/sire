// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "Molecules.pypp.hpp"

namespace bp = boost::python;

#include "SireMol/errors.h"

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "atom.h"

#include "atomselection.h"

#include "chain.h"

#include "cutgroup.h"

#include "editor.hpp"

#include "molecule.h"

#include "molecules.h"

#include "mover.hpp"

#include "partialmolecule.h"

#include "residue.h"

#include "segment.h"

#include "select.h"

#include "selector.hpp"

#include "tostring.h"

#include "viewsofmol.h"

#include <QDebug>

#include <boost/tuple/tuple.hpp>

#include "molecules.h"

SireMol::Molecules __copy__(const SireMol::Molecules &other){ return SireMol::Molecules(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

#include "Helpers/len.hpp"

void register_Molecules_class(){

    { //::SireMol::Molecules
        typedef bp::class_< SireMol::Molecules, bp::bases< SireBase::Property > > Molecules_exposer_t;
        Molecules_exposer_t Molecules_exposer = Molecules_exposer_t( "Molecules", "This class provides a container for lots of molecules. This\nforms a general purpose molecule container, which is used as the argument\nto functions which expect to be passed lots of molecules or parts\nof molecules. This class holds the Molecules using the\nViewsOfMol class, thereby allowing multiple arbitrary views of each\nmolecule to be held.\n\nAuthor: Christopher Woods\n", bp::init< >("Construct an empty set of molecules") );
        bp::scope Molecules_scope( Molecules_exposer );
        Molecules_exposer.def( bp::init< SireMol::MoleculeView const & >(( bp::arg("molecule") ), "Construct a set that contains only the passed view\nof a molecule") );
        Molecules_exposer.def( bp::init< SireMol::ViewsOfMol const & >(( bp::arg("molviews") ), "Construct a set that contains the passed views\nof a molecule") );
        Molecules_exposer.def( bp::init< SireMol::SelectResult const & >(( bp::arg("result") ), "Construct a set that contains the passed search result") );
        Molecules_exposer.def( bp::init< SireMol::Molecules const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMol::Molecules::add
        
            typedef void ( ::SireMol::Molecules::*add_function_type)( ::SireMol::MoleculeView const & ) ;
            add_function_type add_function_value( &::SireMol::Molecules::add );
            
            Molecules_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("molview") )
                , bp::release_gil_policy()
                , "Add the view view to this set of molecules - this\nadds this view even if it already exists in this set" );
        
        }
        { //::SireMol::Molecules::add
        
            typedef void ( ::SireMol::Molecules::*add_function_type)( ::SireMol::ViewsOfMol const & ) ;
            add_function_type add_function_value( &::SireMol::Molecules::add );
            
            Molecules_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("molviews") )
                , bp::release_gil_policy()
                , "Add the views in molviews to this set. This adds\nall of the views, even if they are already present\nin this set." );
        
        }
        { //::SireMol::Molecules::add
        
            typedef void ( ::SireMol::Molecules::*add_function_type)( ::SireMol::Molecules const & ) ;
            add_function_type add_function_value( &::SireMol::Molecules::add );
            
            Molecules_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("molecules") )
                , bp::release_gil_policy()
                , "Add all of the molecules in molecules to this set.\nThis adds all of the molecules, even if they already\nexist in this set" );
        
        }
        { //::SireMol::Molecules::addIfUnique
        
            typedef bool ( ::SireMol::Molecules::*addIfUnique_function_type)( ::SireMol::MoleculeView const & ) ;
            addIfUnique_function_type addIfUnique_function_value( &::SireMol::Molecules::addIfUnique );
            
            Molecules_exposer.def( 
                "addIfUnique"
                , addIfUnique_function_value
                , ( bp::arg("molview") )
                , bp::release_gil_policy()
                , "Add the view view to this set of molecules - this\nadds this view only if it doesnt already exist in\nthis set - this returns whether or not the view\nwas added" );
        
        }
        { //::SireMol::Molecules::addIfUnique
        
            typedef ::SireMol::ViewsOfMol ( ::SireMol::Molecules::*addIfUnique_function_type)( ::SireMol::ViewsOfMol const & ) ;
            addIfUnique_function_type addIfUnique_function_value( &::SireMol::Molecules::addIfUnique );
            
            Molecules_exposer.def( 
                "addIfUnique"
                , addIfUnique_function_value
                , ( bp::arg("molviews") )
                , bp::release_gil_policy()
                , "Add the views in molviews to this set. This only\nadds the views that dont already exist in this set,\nand returns the views that have been added" );
        
        }
        { //::SireMol::Molecules::addIfUnique
        
            typedef ::QList< SireMol::ViewsOfMol > ( ::SireMol::Molecules::*addIfUnique_function_type)( ::SireMol::Molecules const & ) ;
            addIfUnique_function_type addIfUnique_function_value( &::SireMol::Molecules::addIfUnique );
            
            Molecules_exposer.def( 
                "addIfUnique"
                , addIfUnique_function_value
                , ( bp::arg("molecules") )
                , bp::release_gil_policy()
                , "Add all of the views of the molecules in molecules to\nthis set, only if they dont already exist in this set.\nThis returns the views that were added to this set." );
        
        }
        { //::SireMol::Molecules::assertContains
        
            typedef void ( ::SireMol::Molecules::*assertContains_function_type)( ::SireMol::MolNum ) const;
            assertContains_function_type assertContains_function_value( &::SireMol::Molecules::assertContains );
            
            Molecules_exposer.def( 
                "assertContains"
                , assertContains_function_value
                , ( bp::arg("molnum") )
                , bp::release_gil_policy()
                , "Assert that this set contains any of the atoms of\nthe molecule with number molnum\nThrow: SireMol::missing_molecule\n" );
        
        }
        { //::SireMol::Molecules::at
        
            typedef ::SireMol::ViewsOfMol const & ( ::SireMol::Molecules::*at_function_type)( ::SireMol::MolNum ) const;
            at_function_type at_function_value( &::SireMol::Molecules::at );
            
            Molecules_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("molnum") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the view(s) of the molecule that has number molnum\nThrow: SireMol::missing_molecule\n" );
        
        }
        { //::SireMol::Molecules::at
        
            typedef ::SireMol::PartialMolecule ( ::SireMol::Molecules::*at_function_type)( ::boost::tuples::tuple< SireMol::MolNum, SireID::Index, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type > const & ) const;
            at_function_type at_function_value( &::SireMol::Molecules::at );
            
            Molecules_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("viewidx") )
                , bp::release_gil_policy()
                , "Return the view of the molecule identified by molviewidx\nThrow: SireMol::missing_molecule\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMol::Molecules::at
        
            typedef ::SireMol::PartialMolecule ( ::SireMol::Molecules::*at_function_type)( ::SireMol::MolNum,int ) const;
            at_function_type at_function_value( &::SireMol::Molecules::at );
            
            Molecules_exposer.def( 
                "at"
                , at_function_value
                , ( bp::arg("molnum"), bp::arg("idx") )
                , bp::release_gil_policy()
                , "Return the view at index viewidx of the molecule with number\nmolnum from this set\nThrow: SireMol::missing_molecule\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMol::Molecules::back
        
            typedef ::SireMol::ViewsOfMol const & ( ::SireMol::Molecules::*back_function_type)(  ) const;
            back_function_type back_function_value( &::SireMol::Molecules::back );
            
            Molecules_exposer.def( 
                "back"
                , back_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "STL compatible version of Molecules::last()" );
        
        }
        { //::SireMol::Molecules::clear
        
            typedef void ( ::SireMol::Molecules::*clear_function_type)(  ) ;
            clear_function_type clear_function_value( &::SireMol::Molecules::clear );
            
            Molecules_exposer.def( 
                "clear"
                , clear_function_value
                , bp::release_gil_policy()
                , "Completely clear this set of all views of all molecules" );
        
        }
        { //::SireMol::Molecules::contains
        
            typedef bool ( ::SireMol::Molecules::*contains_function_type)( ::SireMol::MolNum ) const;
            contains_function_type contains_function_value( &::SireMol::Molecules::contains );
            
            Molecules_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("molnum") )
                , bp::release_gil_policy()
                , "Return whether or not this set contains any view\nof the molecule at number molnum" );
        
        }
        { //::SireMol::Molecules::contains
        
            typedef bool ( ::SireMol::Molecules::*contains_function_type)( ::SireMol::MoleculeView const & ) const;
            contains_function_type contains_function_value( &::SireMol::Molecules::contains );
            
            Molecules_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("molview") )
                , bp::release_gil_policy()
                , "Return whether or not this set contains the view molview\n- note that this means that it actually contains this\nspecific view" );
        
        }
        { //::SireMol::Molecules::contains
        
            typedef bool ( ::SireMol::Molecules::*contains_function_type)( ::SireMol::ViewsOfMol const & ) const;
            contains_function_type contains_function_value( &::SireMol::Molecules::contains );
            
            Molecules_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("molviews") )
                , bp::release_gil_policy()
                , "Return whether or not this set contains all of the views\nin molviews - note that this means that it actually\ncontains each specific view" );
        
        }
        { //::SireMol::Molecules::contains
        
            typedef bool ( ::SireMol::Molecules::*contains_function_type)( ::SireMol::Molecules const & ) const;
            contains_function_type contains_function_value( &::SireMol::Molecules::contains );
            
            Molecules_exposer.def( 
                "contains"
                , contains_function_value
                , ( bp::arg("molecules") )
                , bp::release_gil_policy()
                , "Return whether or not this set contains all of the views\nof all of the molecules in molecules" );
        
        }
        { //::SireMol::Molecules::count
        
            typedef int ( ::SireMol::Molecules::*count_function_type)(  ) const;
            count_function_type count_function_value( &::SireMol::Molecules::count );
            
            Molecules_exposer.def( 
                "count"
                , count_function_value
                , bp::release_gil_policy()
                , "Return the number of molecules in this set" );
        
        }
        { //::SireMol::Molecules::deleteFrame
        
            typedef void ( ::SireMol::Molecules::*deleteFrame_function_type)( int ) ;
            deleteFrame_function_type deleteFrame_function_value( &::SireMol::Molecules::deleteFrame );
            
            Molecules_exposer.def( 
                "deleteFrame"
                , deleteFrame_function_value
                , ( bp::arg("frame") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Molecules::deleteFrame
        
            typedef void ( ::SireMol::Molecules::*deleteFrame_function_type)( int,::SireBase::PropertyMap const & ) ;
            deleteFrame_function_type deleteFrame_function_value( &::SireMol::Molecules::deleteFrame );
            
            Molecules_exposer.def( 
                "deleteFrame"
                , deleteFrame_function_value
                , ( bp::arg("frame"), bp::arg("map") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Molecules::first
        
            typedef ::SireMol::ViewsOfMol const & ( ::SireMol::Molecules::*first_function_type)(  ) const;
            first_function_type first_function_value( &::SireMol::Molecules::first );
            
            Molecules_exposer.def( 
                "first"
                , first_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return a reference to the first molecule in this set.\nThis throws an exception if this set is empty.\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMol::Molecules::front
        
            typedef ::SireMol::ViewsOfMol const & ( ::SireMol::Molecules::*front_function_type)(  ) const;
            front_function_type front_function_value( &::SireMol::Molecules::front );
            
            Molecules_exposer.def( 
                "front"
                , front_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "STL compatible version of Molecules::first()" );
        
        }
        { //::SireMol::Molecules::intersects
        
            typedef bool ( ::SireMol::Molecules::*intersects_function_type)( ::SireMol::MoleculeView const & ) const;
            intersects_function_type intersects_function_value( &::SireMol::Molecules::intersects );
            
            Molecules_exposer.def( 
                "intersects"
                , intersects_function_value
                , ( bp::arg("molview") )
                , bp::release_gil_policy()
                , "Return whether or not this set intersects with the view molview\n- this returns true if any of the atoms in molview are\nalso present in any of the views of that molecule in this set" );
        
        }
        { //::SireMol::Molecules::intersects
        
            typedef bool ( ::SireMol::Molecules::*intersects_function_type)( ::SireMol::Molecules const & ) const;
            intersects_function_type intersects_function_value( &::SireMol::Molecules::intersects );
            
            Molecules_exposer.def( 
                "intersects"
                , intersects_function_value
                , ( bp::arg("other") )
                , bp::release_gil_policy()
                , "Return whether or not this set intersects with molecules -\nthis returns true if any of the atoms in any of the views\nof any of the molecules in molecules are also present in any\nof the views of any of the molecules in this set" );
        
        }
        { //::SireMol::Molecules::isEmpty
        
            typedef bool ( ::SireMol::Molecules::*isEmpty_function_type)(  ) const;
            isEmpty_function_type isEmpty_function_value( &::SireMol::Molecules::isEmpty );
            
            Molecules_exposer.def( 
                "isEmpty"
                , isEmpty_function_value
                , bp::release_gil_policy()
                , "Return whether this set is empty" );
        
        }
        { //::SireMol::Molecules::last
        
            typedef ::SireMol::ViewsOfMol const & ( ::SireMol::Molecules::*last_function_type)(  ) const;
            last_function_type last_function_value( &::SireMol::Molecules::last );
            
            Molecules_exposer.def( 
                "last"
                , last_function_value
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return a reference to the last molecule in this set.\nThis throws an exception if this set is empty.\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMol::Molecules::loadFrame
        
            typedef void ( ::SireMol::Molecules::*loadFrame_function_type)( int ) ;
            loadFrame_function_type loadFrame_function_value( &::SireMol::Molecules::loadFrame );
            
            Molecules_exposer.def( 
                "loadFrame"
                , loadFrame_function_value
                , ( bp::arg("frame") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Molecules::loadFrame
        
            typedef void ( ::SireMol::Molecules::*loadFrame_function_type)( int,::SireBase::PropertyMap const & ) ;
            loadFrame_function_type loadFrame_function_value( &::SireMol::Molecules::loadFrame );
            
            Molecules_exposer.def( 
                "loadFrame"
                , loadFrame_function_value
                , ( bp::arg("frame"), bp::arg("map") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Molecules::molNums
        
            typedef ::QSet< SireMol::MolNum > ( ::SireMol::Molecules::*molNums_function_type)(  ) const;
            molNums_function_type molNums_function_value( &::SireMol::Molecules::molNums );
            
            Molecules_exposer.def( 
                "molNums"
                , molNums_function_value
                , bp::release_gil_policy()
                , "Return the numbers of all of the molecules in this set" );
        
        }
        { //::SireMol::Molecules::molecule
        
            typedef ::SireMol::ViewsOfMol const & ( ::SireMol::Molecules::*molecule_function_type)( ::SireMol::MolNum ) const;
            molecule_function_type molecule_function_value( &::SireMol::Molecules::molecule );
            
            Molecules_exposer.def( 
                "molecule"
                , molecule_function_value
                , ( bp::arg("molnum") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "Return the molecule at number molnum\nThrow: SireMol::missing_molecule\n" );
        
        }
        { //::SireMol::Molecules::nFrames
        
            typedef int ( ::SireMol::Molecules::*nFrames_function_type)(  ) const;
            nFrames_function_type nFrames_function_value( &::SireMol::Molecules::nFrames );
            
            Molecules_exposer.def( 
                "nFrames"
                , nFrames_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Molecules::nFrames
        
            typedef int ( ::SireMol::Molecules::*nFrames_function_type)( ::SireBase::PropertyMap const & ) const;
            nFrames_function_type nFrames_function_value( &::SireMol::Molecules::nFrames );
            
            Molecules_exposer.def( 
                "nFrames"
                , nFrames_function_value
                , ( bp::arg("map") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Molecules::nMolecules
        
            typedef int ( ::SireMol::Molecules::*nMolecules_function_type)(  ) const;
            nMolecules_function_type nMolecules_function_value( &::SireMol::Molecules::nMolecules );
            
            Molecules_exposer.def( 
                "nMolecules"
                , nMolecules_function_value
                , bp::release_gil_policy()
                , "Return the number of molecules in this set" );
        
        }
        { //::SireMol::Molecules::nViews
        
            typedef int ( ::SireMol::Molecules::*nViews_function_type)(  ) const;
            nViews_function_type nViews_function_value( &::SireMol::Molecules::nViews );
            
            Molecules_exposer.def( 
                "nViews"
                , nViews_function_value
                , bp::release_gil_policy()
                , "Return the total number of views in this set" );
        
        }
        { //::SireMol::Molecules::nViews
        
            typedef int ( ::SireMol::Molecules::*nViews_function_type)( ::SireMol::MolNum ) const;
            nViews_function_type nViews_function_value( &::SireMol::Molecules::nViews );
            
            Molecules_exposer.def( 
                "nViews"
                , nViews_function_value
                , ( bp::arg("molnum") )
                , bp::release_gil_policy()
                , "Return the number of views of the molecule with number molnum\nThrow: SireMol::missing_molecule\n" );
        
        }
        Molecules_exposer.def( bp::self != bp::self );
        Molecules_exposer.def( bp::self + bp::self );
        Molecules_exposer.def( bp::self - bp::self );
        { //::SireMol::Molecules::operator=
        
            typedef ::SireMol::Molecules & ( ::SireMol::Molecules::*assign_function_type)( ::SireMol::Molecules const & ) ;
            assign_function_type assign_function_value( &::SireMol::Molecules::operator= );
            
            Molecules_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        Molecules_exposer.def( bp::self == bp::self );
        { //::SireMol::Molecules::operator[]
        
            typedef ::SireMol::ViewsOfMol const & ( ::SireMol::Molecules::*__getitem___function_type)( ::SireMol::MolNum ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::Molecules::operator[] );
            
            Molecules_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("molnum") )
                , bp::return_value_policy< bp::copy_const_reference >()
                , "" );
        
        }
        { //::SireMol::Molecules::operator[]
        
            typedef ::SireMol::PartialMolecule ( ::SireMol::Molecules::*__getitem___function_type)( ::boost::tuples::tuple< SireMol::MolNum, SireID::Index, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type, boost::tuples::null_type > const & ) const;
            __getitem___function_type __getitem___function_value( &::SireMol::Molecules::operator[] );
            
            Molecules_exposer.def( 
                "__getitem__"
                , __getitem___function_value
                , ( bp::arg("viewidx") )
                , "" );
        
        }
        { //::SireMol::Molecules::remove
        
            typedef bool ( ::SireMol::Molecules::*remove_function_type)( ::SireMol::MoleculeView const & ) ;
            remove_function_type remove_function_value( &::SireMol::Molecules::remove );
            
            Molecules_exposer.def( 
                "remove"
                , remove_function_value
                , ( bp::arg("molview") )
                , bp::release_gil_policy()
                , "Remove the view molview from this set. This only\nremoves the first copy of this view from this set\n(if multiple copies are present), and returns\nwhether or not a view was removed" );
        
        }
        { //::SireMol::Molecules::remove
        
            typedef ::SireMol::ViewsOfMol ( ::SireMol::Molecules::*remove_function_type)( ::SireMol::ViewsOfMol const & ) ;
            remove_function_type remove_function_value( &::SireMol::Molecules::remove );
            
            Molecules_exposer.def( 
                "remove"
                , remove_function_value
                , ( bp::arg("molviews") )
                , bp::release_gil_policy()
                , "Remove all of the views in molviews from this set. This\nonly removes the first copy of the view if multiple copies\nexist, and returns the views that were successfully removed." );
        
        }
        { //::SireMol::Molecules::remove
        
            typedef ::QList< SireMol::ViewsOfMol > ( ::SireMol::Molecules::*remove_function_type)( ::SireMol::Molecules const & ) ;
            remove_function_type remove_function_value( &::SireMol::Molecules::remove );
            
            Molecules_exposer.def( 
                "remove"
                , remove_function_value
                , ( bp::arg("molecules") )
                , bp::release_gil_policy()
                , "Remove all of the views of all of the molecules in molecules.\nThis only removes the first copy of any views that appear\nmultiple times in this set. This returns all of the views that\nwere successfully removed" );
        
        }
        { //::SireMol::Molecules::remove
        
            typedef ::SireMol::ViewsOfMol ( ::SireMol::Molecules::*remove_function_type)( ::SireMol::MolNum ) ;
            remove_function_type remove_function_value( &::SireMol::Molecules::remove );
            
            Molecules_exposer.def( 
                "remove"
                , remove_function_value
                , ( bp::arg("molnum") )
                , bp::release_gil_policy()
                , "Remove all views of the molecule with number molnum. This\nreturns the views of the molecule in this set, if it exists,\nor an empty set of views if it doesnt exist." );
        
        }
        { //::SireMol::Molecules::removeAll
        
            typedef bool ( ::SireMol::Molecules::*removeAll_function_type)( ::SireMol::MoleculeView const & ) ;
            removeAll_function_type removeAll_function_value( &::SireMol::Molecules::removeAll );
            
            Molecules_exposer.def( 
                "removeAll"
                , removeAll_function_value
                , ( bp::arg("molview") )
                , bp::release_gil_policy()
                , "Remove all copies of the view view from this set.\nThis returns whether any views were removed" );
        
        }
        { //::SireMol::Molecules::removeAll
        
            typedef ::SireMol::ViewsOfMol ( ::SireMol::Molecules::*removeAll_function_type)( ::SireMol::ViewsOfMol const & ) ;
            removeAll_function_type removeAll_function_value( &::SireMol::Molecules::removeAll );
            
            Molecules_exposer.def( 
                "removeAll"
                , removeAll_function_value
                , ( bp::arg("molviews") )
                , bp::release_gil_policy()
                , "Remove all of the views in molviews from this set. This\nremoves all copies of the views if multiple copies exist,\nand returns the views that were successfully removed." );
        
        }
        { //::SireMol::Molecules::removeAll
        
            typedef ::QList< SireMol::ViewsOfMol > ( ::SireMol::Molecules::*removeAll_function_type)( ::SireMol::Molecules const & ) ;
            removeAll_function_type removeAll_function_value( &::SireMol::Molecules::removeAll );
            
            Molecules_exposer.def( 
                "removeAll"
                , removeAll_function_value
                , ( bp::arg("molecules") )
                , bp::release_gil_policy()
                , "Remove all of the views of all of the molecules in molecules.\nThis removes all copies of any views that appear\nmultiple times in this set. This returns all of the views that\nwere successfully removed" );
        
        }
        { //::SireMol::Molecules::removeAll
        
            typedef bool ( ::SireMol::Molecules::*removeAll_function_type)(  ) ;
            removeAll_function_type removeAll_function_value( &::SireMol::Molecules::removeAll );
            
            Molecules_exposer.def( 
                "removeAll"
                , removeAll_function_value
                , bp::release_gil_policy()
                , "Remove all views of all molecules from this set. This returns\nwhether or not this changes this set" );
        
        }
        { //::SireMol::Molecules::removeDuplicates
        
            typedef bool ( ::SireMol::Molecules::*removeDuplicates_function_type)(  ) ;
            removeDuplicates_function_type removeDuplicates_function_value( &::SireMol::Molecules::removeDuplicates );
            
            Molecules_exposer.def( 
                "removeDuplicates"
                , removeDuplicates_function_value
                , bp::release_gil_policy()
                , "This removes all duplicated views from this set. This returns\nwhether or not this changes the set (whether or not there\nwere any duplicates)" );
        
        }
        { //::SireMol::Molecules::reserve
        
            typedef void ( ::SireMol::Molecules::*reserve_function_type)( int ) ;
            reserve_function_type reserve_function_value( &::SireMol::Molecules::reserve );
            
            Molecules_exposer.def( 
                "reserve"
                , reserve_function_value
                , ( bp::arg("nmolecules") )
                , bp::release_gil_policy()
                , "Reserve enough space for nmolecules molecules. This\nwill reserve the memory so that reallocations are minimised" );
        
        }
        { //::SireMol::Molecules::saveFrame
        
            typedef void ( ::SireMol::Molecules::*saveFrame_function_type)( int ) ;
            saveFrame_function_type saveFrame_function_value( &::SireMol::Molecules::saveFrame );
            
            Molecules_exposer.def( 
                "saveFrame"
                , saveFrame_function_value
                , ( bp::arg("frame") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Molecules::saveFrame
        
            typedef void ( ::SireMol::Molecules::*saveFrame_function_type)(  ) ;
            saveFrame_function_type saveFrame_function_value( &::SireMol::Molecules::saveFrame );
            
            Molecules_exposer.def( 
                "saveFrame"
                , saveFrame_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Molecules::saveFrame
        
            typedef void ( ::SireMol::Molecules::*saveFrame_function_type)( int,::SireBase::PropertyMap const & ) ;
            saveFrame_function_type saveFrame_function_value( &::SireMol::Molecules::saveFrame );
            
            Molecules_exposer.def( 
                "saveFrame"
                , saveFrame_function_value
                , ( bp::arg("frame"), bp::arg("map") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Molecules::saveFrame
        
            typedef void ( ::SireMol::Molecules::*saveFrame_function_type)( ::SireBase::PropertyMap const & ) ;
            saveFrame_function_type saveFrame_function_value( &::SireMol::Molecules::saveFrame );
            
            Molecules_exposer.def( 
                "saveFrame"
                , saveFrame_function_value
                , ( bp::arg("map") )
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Molecules::search
        
            typedef ::SireMol::SelectResult ( ::SireMol::Molecules::*search_function_type)( ::QString const & ) const;
            search_function_type search_function_value( &::SireMol::Molecules::search );
            
            Molecules_exposer.def( 
                "search"
                , search_function_value
                , ( bp::arg("search_string") )
                , bp::release_gil_policy()
                , "Return the result of searching these molecules with search_string" );
        
        }
        { //::SireMol::Molecules::toString
        
            typedef ::QString ( ::SireMol::Molecules::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::Molecules::toString );
            
            Molecules_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "Return a string representation of this set of molecules" );
        
        }
        { //::SireMol::Molecules::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::Molecules::typeName );
            
            Molecules_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::Molecules::unite
        
            typedef bool ( ::SireMol::Molecules::*unite_function_type)( ::SireMol::MoleculeView const & ) ;
            unite_function_type unite_function_value( &::SireMol::Molecules::unite );
            
            Molecules_exposer.def( 
                "unite"
                , unite_function_value
                , ( bp::arg("molview") )
                , bp::release_gil_policy()
                , "Synonym for Molecules::addIfUnique(molview)" );
        
        }
        { //::SireMol::Molecules::unite
        
            typedef ::SireMol::ViewsOfMol ( ::SireMol::Molecules::*unite_function_type)( ::SireMol::ViewsOfMol const & ) ;
            unite_function_type unite_function_value( &::SireMol::Molecules::unite );
            
            Molecules_exposer.def( 
                "unite"
                , unite_function_value
                , ( bp::arg("molviews") )
                , bp::release_gil_policy()
                , "Synonym for Molecules::addIfUnique(molviews)" );
        
        }
        { //::SireMol::Molecules::unite
        
            typedef ::QList< SireMol::ViewsOfMol > ( ::SireMol::Molecules::*unite_function_type)( ::SireMol::Molecules const & ) ;
            unite_function_type unite_function_value( &::SireMol::Molecules::unite );
            
            Molecules_exposer.def( 
                "unite"
                , unite_function_value
                , ( bp::arg("other") )
                , bp::release_gil_policy()
                , "Synonym for Molecules::addIfUnique(molecules)" );
        
        }
        { //::SireMol::Molecules::uniteViews
        
            typedef bool ( ::SireMol::Molecules::*uniteViews_function_type)(  ) ;
            uniteViews_function_type uniteViews_function_value( &::SireMol::Molecules::uniteViews );
            
            Molecules_exposer.def( 
                "uniteViews"
                , uniteViews_function_value
                , bp::release_gil_policy()
                , "Unite all of the views in this set so that each molecule has\nonly a single view that is the union of all of its views.\nReturn whether or not this changes the set." );
        
        }
        { //::SireMol::Molecules::update
        
            typedef bool ( ::SireMol::Molecules::*update_function_type)( ::SireMol::MoleculeData const & ) ;
            update_function_type update_function_value( &::SireMol::Molecules::update );
            
            Molecules_exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("moldata") )
                , bp::release_gil_policy()
                , "Update the views of the molecule whose data is in moldata\nin this set so that it is at the same version as moldata.\nThis does nothing if this molecule is not in this set.\nThis returns whether or not the molecule was updated" );
        
        }
        { //::SireMol::Molecules::update
        
            typedef bool ( ::SireMol::Molecules::*update_function_type)( ::SireMol::MoleculeView const & ) ;
            update_function_type update_function_value( &::SireMol::Molecules::update );
            
            Molecules_exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("molview") )
                , bp::release_gil_policy()
                , "Update the views of the molecule viewed by molview\nin this set so that they have the same molecule version\nas molview. This returns whether or not the molecule\nwas updated." );
        
        }
        { //::SireMol::Molecules::update
        
            typedef ::QList< SireMol::Molecule > ( ::SireMol::Molecules::*update_function_type)( ::SireMol::Molecules const & ) ;
            update_function_type update_function_value( &::SireMol::Molecules::update );
            
            Molecules_exposer.def( 
                "update"
                , update_function_value
                , ( bp::arg("molecules") )
                , bp::release_gil_policy()
                , "Update the views in this set so that they have the\nsame molecule versions as the molecules in molecules.\nThis returns the molecules that have been updated." );
        
        }
        { //::SireMol::Molecules::what
        
            typedef char const * ( ::SireMol::Molecules::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMol::Molecules::what );
            
            Molecules_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        Molecules_exposer.staticmethod( "typeName" );
        Molecules_exposer.def( "__copy__", &__copy__);
        Molecules_exposer.def( "__deepcopy__", &__copy__);
        Molecules_exposer.def( "clone", &__copy__);
        Molecules_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::Molecules >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Molecules_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::Molecules >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        Molecules_exposer.def_pickle(sire_pickle_suite< ::SireMol::Molecules >());
        Molecules_exposer.def( "__str__", &__str__< ::SireMol::Molecules > );
        Molecules_exposer.def( "__repr__", &__str__< ::SireMol::Molecules > );
        Molecules_exposer.def( "__len__", &__len_count< ::SireMol::Molecules > );
    }

}
