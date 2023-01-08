// This file has been generated by Py++.

// (C) Christopher Woods, GPL >= 3 License

#include "boost/python.hpp"
#include "Helpers/clone_const_reference.hpp"
#include "MolStructureEditor.pypp.hpp"

namespace bp = boost::python;

#include "SireStream/datastream.h"

#include "SireStream/shareddatastream.h"

#include "atom.h"

#include "atomeditor.h"

#include "cgeditor.h"

#include "chain.h"

#include "chaineditor.h"

#include "cutgroup.h"

#include "moleditor.h"

#include "mover.hpp"

#include "reseditor.h"

#include "residue.h"

#include "segeditor.h"

#include "segment.h"

#include "selector.hpp"

#include "moleditor.h"

SireMol::MolStructureEditor __copy__(const SireMol::MolStructureEditor &other){ return SireMol::MolStructureEditor(other); }

#include "Qt/qdatastream.hpp"

#include "Helpers/str.hpp"

#include "Helpers/release_gil_policy.hpp"

void register_MolStructureEditor_class(){

    { //::SireMol::MolStructureEditor
        typedef bp::class_< SireMol::MolStructureEditor > MolStructureEditor_exposer_t;
        MolStructureEditor_exposer_t MolStructureEditor_exposer = MolStructureEditor_exposer_t( "MolStructureEditor", "This class is used to edit structural parts of the molecule,\ni.e. adding, moving or removing atoms, residues etc.\n\nAuthor: Christopher Woods\n", bp::init< >("Null constructor") );
        bp::scope MolStructureEditor_scope( MolStructureEditor_exposer );
        MolStructureEditor_exposer.def( bp::init< SireMol::MoleculeView const & >(( bp::arg("molview") ), "Construct a structure editor of the molecule viewed by molview") );
        MolStructureEditor_exposer.def( bp::init< SireMol::StructureEditor const & >(( bp::arg("other") ), "Copy constructor") );
        MolStructureEditor_exposer.def( bp::init< SireMol::MolStructureEditor const & >(( bp::arg("other") ), "Copy constructor") );
        { //::SireMol::MolStructureEditor::add
        
            typedef ::SireMol::AtomStructureEditor ( ::SireMol::MolStructureEditor::*add_function_type)( ::SireMol::AtomName const & ) ;
            add_function_type add_function_value( &::SireMol::MolStructureEditor::add );
            
            MolStructureEditor_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("atom") )
                , bp::release_gil_policy()
                , "Add an atom called name to this molecule and return an\neditor for this atom" );
        
        }
        { //::SireMol::MolStructureEditor::add
        
            typedef ::SireMol::AtomStructureEditor ( ::SireMol::MolStructureEditor::*add_function_type)( ::SireMol::AtomNum const & ) ;
            add_function_type add_function_value( &::SireMol::MolStructureEditor::add );
            
            MolStructureEditor_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("atom") )
                , bp::release_gil_policy()
                , "Add an atom called name to this molecule and return an\neditor for this atom" );
        
        }
        { //::SireMol::MolStructureEditor::add
        
            typedef ::SireMol::ResStructureEditor ( ::SireMol::MolStructureEditor::*add_function_type)( ::SireMol::ResName const & ) ;
            add_function_type add_function_value( &::SireMol::MolStructureEditor::add );
            
            MolStructureEditor_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("residue") )
                , bp::release_gil_policy()
                , "Add a residue called name to this molecule and return an\neditor for this residue" );
        
        }
        { //::SireMol::MolStructureEditor::add
        
            typedef ::SireMol::ResStructureEditor ( ::SireMol::MolStructureEditor::*add_function_type)( ::SireMol::ResNum const & ) ;
            add_function_type add_function_value( &::SireMol::MolStructureEditor::add );
            
            MolStructureEditor_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("residue") )
                , bp::release_gil_policy()
                , "Add a residue with number number to this molecule and\nreturn an editor for this residue" );
        
        }
        { //::SireMol::MolStructureEditor::add
        
            typedef ::SireMol::CGStructureEditor ( ::SireMol::MolStructureEditor::*add_function_type)( ::SireMol::CGName const & ) ;
            add_function_type add_function_value( &::SireMol::MolStructureEditor::add );
            
            MolStructureEditor_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("cutgroup") )
                , bp::release_gil_policy()
                , "Add a CutGroup called name to this molecule and return an\neditor for this CutGroup" );
        
        }
        { //::SireMol::MolStructureEditor::add
        
            typedef ::SireMol::ChainStructureEditor ( ::SireMol::MolStructureEditor::*add_function_type)( ::SireMol::ChainName const & ) ;
            add_function_type add_function_value( &::SireMol::MolStructureEditor::add );
            
            MolStructureEditor_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("chain") )
                , bp::release_gil_policy()
                , "Add a chain called name to this molecule and return an\neditor for this chain" );
        
        }
        { //::SireMol::MolStructureEditor::add
        
            typedef ::SireMol::SegStructureEditor ( ::SireMol::MolStructureEditor::*add_function_type)( ::SireMol::SegName const & ) ;
            add_function_type add_function_value( &::SireMol::MolStructureEditor::add );
            
            MolStructureEditor_exposer.def( 
                "add"
                , add_function_value
                , ( bp::arg("segment") )
                , bp::release_gil_policy()
                , "Add a segment called name to this molecule and return an\neditor for this segment" );
        
        }
        { //::SireMol::MolStructureEditor::atom
        
            typedef ::SireMol::AtomStructureEditor ( ::SireMol::MolStructureEditor::*atom_function_type)( ::SireMol::AtomID const & ) ;
            atom_function_type atom_function_value( &::SireMol::MolStructureEditor::atom );
            
            MolStructureEditor_exposer.def( 
                "atom"
                , atom_function_value
                , ( bp::arg("atomid") )
                , bp::release_gil_policy()
                , "Return an editor for the atom at ID atomid\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMol::MolStructureEditor::chain
        
            typedef ::SireMol::ChainStructureEditor ( ::SireMol::MolStructureEditor::*chain_function_type)( ::SireMol::ChainID const & ) ;
            chain_function_type chain_function_value( &::SireMol::MolStructureEditor::chain );
            
            MolStructureEditor_exposer.def( 
                "chain"
                , chain_function_value
                , ( bp::arg("chainid") )
                , bp::release_gil_policy()
                , "Return an editor for the chain at ID chainid\nThrow: SireMol::missing_chain\nThrow: SireMol::duplicate_chain\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMol::MolStructureEditor::commit
        
            typedef ::SireMol::Molecule ( ::SireMol::MolStructureEditor::*commit_function_type)(  ) const;
            commit_function_type commit_function_value( &::SireMol::MolStructureEditor::commit );
            
            MolStructureEditor_exposer.def( 
                "commit"
                , commit_function_value
                , bp::release_gil_policy()
                , "Commit the changes and return a Molecule copy" );
        
        }
        { //::SireMol::MolStructureEditor::cutGroup
        
            typedef ::SireMol::CGStructureEditor ( ::SireMol::MolStructureEditor::*cutGroup_function_type)( ::SireMol::CGID const & ) ;
            cutGroup_function_type cutGroup_function_value( &::SireMol::MolStructureEditor::cutGroup );
            
            MolStructureEditor_exposer.def( 
                "cutGroup"
                , cutGroup_function_value
                , ( bp::arg("cgid") )
                , bp::release_gil_policy()
                , "Return an editor for the CutGroup at ID cgid\nThrow: SireMol::missing_cutgroup\nThrow: SireMol::duplicate_cutgroup\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMol::MolStructureEditor::nAtoms
        
            typedef int ( ::SireMol::MolStructureEditor::*nAtoms_function_type)(  ) const;
            nAtoms_function_type nAtoms_function_value( &::SireMol::MolStructureEditor::nAtoms );
            
            MolStructureEditor_exposer.def( 
                "nAtoms"
                , nAtoms_function_value
                , bp::release_gil_policy()
                , "Return the number of atoms in this molecule (may be zero)" );
        
        }
        { //::SireMol::MolStructureEditor::nChains
        
            typedef int ( ::SireMol::MolStructureEditor::*nChains_function_type)(  ) const;
            nChains_function_type nChains_function_value( &::SireMol::MolStructureEditor::nChains );
            
            MolStructureEditor_exposer.def( 
                "nChains"
                , nChains_function_value
                , bp::release_gil_policy()
                , "Return the number of chains in this molecule (may be zero)" );
        
        }
        { //::SireMol::MolStructureEditor::nCutGroups
        
            typedef int ( ::SireMol::MolStructureEditor::*nCutGroups_function_type)(  ) const;
            nCutGroups_function_type nCutGroups_function_value( &::SireMol::MolStructureEditor::nCutGroups );
            
            MolStructureEditor_exposer.def( 
                "nCutGroups"
                , nCutGroups_function_value
                , bp::release_gil_policy()
                , "Return the number of CutGroups in this molecule (may be zero)" );
        
        }
        { //::SireMol::MolStructureEditor::nResidues
        
            typedef int ( ::SireMol::MolStructureEditor::*nResidues_function_type)(  ) const;
            nResidues_function_type nResidues_function_value( &::SireMol::MolStructureEditor::nResidues );
            
            MolStructureEditor_exposer.def( 
                "nResidues"
                , nResidues_function_value
                , bp::release_gil_policy()
                , "Return the number of residues in this molecule (may be zero)" );
        
        }
        { //::SireMol::MolStructureEditor::nSegments
        
            typedef int ( ::SireMol::MolStructureEditor::*nSegments_function_type)(  ) const;
            nSegments_function_type nSegments_function_value( &::SireMol::MolStructureEditor::nSegments );
            
            MolStructureEditor_exposer.def( 
                "nSegments"
                , nSegments_function_value
                , bp::release_gil_policy()
                , "Return the number of segments in this molecule (may be zero)" );
        
        }
        { //::SireMol::MolStructureEditor::name
        
            typedef ::SireMol::MolName const & ( ::SireMol::MolStructureEditor::*name_function_type)(  ) const;
            name_function_type name_function_value( &::SireMol::MolStructureEditor::name );
            
            MolStructureEditor_exposer.def( 
                "name"
                , name_function_value
                , bp::return_value_policy<bp::clone_const_reference, bp::release_gil_policy>()
                , "Return the name of this molecule" );
        
        }
        { //::SireMol::MolStructureEditor::number
        
            typedef ::SireMol::MolNum ( ::SireMol::MolStructureEditor::*number_function_type)(  ) const;
            number_function_type number_function_value( &::SireMol::MolStructureEditor::number );
            
            MolStructureEditor_exposer.def( 
                "number"
                , number_function_value
                , bp::release_gil_policy()
                , "Return this molecules ID number" );
        
        }
        { //::SireMol::MolStructureEditor::operator=
        
            typedef ::SireMol::MolStructureEditor & ( ::SireMol::MolStructureEditor::*assign_function_type)( ::SireMol::MoleculeView const & ) ;
            assign_function_type assign_function_value( &::SireMol::MolStructureEditor::operator= );
            
            MolStructureEditor_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("molview") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::MolStructureEditor::operator=
        
            typedef ::SireMol::MolStructureEditor & ( ::SireMol::MolStructureEditor::*assign_function_type)( ::SireMol::StructureEditor const & ) ;
            assign_function_type assign_function_value( &::SireMol::MolStructureEditor::operator= );
            
            MolStructureEditor_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::MolStructureEditor::operator=
        
            typedef ::SireMol::MolStructureEditor & ( ::SireMol::MolStructureEditor::*assign_function_type)( ::SireMol::MolStructureEditor const & ) ;
            assign_function_type assign_function_value( &::SireMol::MolStructureEditor::operator= );
            
            MolStructureEditor_exposer.def( 
                "assign"
                , assign_function_value
                , ( bp::arg("other") )
                , bp::return_self< >()
                , "" );
        
        }
        { //::SireMol::MolStructureEditor::remove
        
            typedef ::SireMol::MolStructureEditor & ( ::SireMol::MolStructureEditor::*remove_function_type)( ::SireMol::AtomID const & ) ;
            remove_function_type remove_function_value( &::SireMol::MolStructureEditor::remove );
            
            MolStructureEditor_exposer.def( 
                "remove"
                , remove_function_value
                , ( bp::arg("atomid") )
                , bp::return_self< >()
                , "Remove all atoms that match the ID atomid - this does\nnothing if there are no atoms that match this ID" );
        
        }
        { //::SireMol::MolStructureEditor::remove
        
            typedef ::SireMol::MolStructureEditor & ( ::SireMol::MolStructureEditor::*remove_function_type)( ::SireMol::CGID const & ) ;
            remove_function_type remove_function_value( &::SireMol::MolStructureEditor::remove );
            
            MolStructureEditor_exposer.def( 
                "remove"
                , remove_function_value
                , ( bp::arg("cgid") )
                , bp::return_self< >()
                , "Remove all CutGroups that match the ID cgid - this does\nnothing if there are no CutGroups that match this ID" );
        
        }
        { //::SireMol::MolStructureEditor::remove
        
            typedef ::SireMol::MolStructureEditor & ( ::SireMol::MolStructureEditor::*remove_function_type)( ::SireMol::ResID const & ) ;
            remove_function_type remove_function_value( &::SireMol::MolStructureEditor::remove );
            
            MolStructureEditor_exposer.def( 
                "remove"
                , remove_function_value
                , ( bp::arg("resid") )
                , bp::return_self< >()
                , "Remove all residues that match the ID resid - this does\nnothing if there are no residues that match this ID" );
        
        }
        { //::SireMol::MolStructureEditor::remove
        
            typedef ::SireMol::MolStructureEditor & ( ::SireMol::MolStructureEditor::*remove_function_type)( ::SireMol::ChainID const & ) ;
            remove_function_type remove_function_value( &::SireMol::MolStructureEditor::remove );
            
            MolStructureEditor_exposer.def( 
                "remove"
                , remove_function_value
                , ( bp::arg("chainid") )
                , bp::return_self< >()
                , "Remove all chains that match the ID chainid - this does\nnothing if there are no chains that match this ID" );
        
        }
        { //::SireMol::MolStructureEditor::remove
        
            typedef ::SireMol::MolStructureEditor & ( ::SireMol::MolStructureEditor::*remove_function_type)( ::SireMol::SegID const & ) ;
            remove_function_type remove_function_value( &::SireMol::MolStructureEditor::remove );
            
            MolStructureEditor_exposer.def( 
                "remove"
                , remove_function_value
                , ( bp::arg("segid") )
                , bp::return_self< >()
                , "Remove all segments that match the ID segid - this does\nnothing if there are no segments that match this ID" );
        
        }
        { //::SireMol::MolStructureEditor::removeAllAtoms
        
            typedef ::SireMol::MolStructureEditor & ( ::SireMol::MolStructureEditor::*removeAllAtoms_function_type)(  ) ;
            removeAllAtoms_function_type removeAllAtoms_function_value( &::SireMol::MolStructureEditor::removeAllAtoms );
            
            MolStructureEditor_exposer.def( 
                "removeAllAtoms"
                , removeAllAtoms_function_value
                , bp::return_self< >()
                , "Remove all atoms from this molecule" );
        
        }
        { //::SireMol::MolStructureEditor::removeAllChains
        
            typedef ::SireMol::MolStructureEditor & ( ::SireMol::MolStructureEditor::*removeAllChains_function_type)(  ) ;
            removeAllChains_function_type removeAllChains_function_value( &::SireMol::MolStructureEditor::removeAllChains );
            
            MolStructureEditor_exposer.def( 
                "removeAllChains"
                , removeAllChains_function_value
                , bp::return_self< >()
                , "Remove all chains from this molecule" );
        
        }
        { //::SireMol::MolStructureEditor::removeAllCutGroups
        
            typedef ::SireMol::MolStructureEditor & ( ::SireMol::MolStructureEditor::*removeAllCutGroups_function_type)(  ) ;
            removeAllCutGroups_function_type removeAllCutGroups_function_value( &::SireMol::MolStructureEditor::removeAllCutGroups );
            
            MolStructureEditor_exposer.def( 
                "removeAllCutGroups"
                , removeAllCutGroups_function_value
                , bp::return_self< >()
                , "Remove all CutGroups from this molecule" );
        
        }
        { //::SireMol::MolStructureEditor::removeAllResidues
        
            typedef ::SireMol::MolStructureEditor & ( ::SireMol::MolStructureEditor::*removeAllResidues_function_type)(  ) ;
            removeAllResidues_function_type removeAllResidues_function_value( &::SireMol::MolStructureEditor::removeAllResidues );
            
            MolStructureEditor_exposer.def( 
                "removeAllResidues"
                , removeAllResidues_function_value
                , bp::return_self< >()
                , "Remove all residues from this molecule" );
        
        }
        { //::SireMol::MolStructureEditor::removeAllSegments
        
            typedef ::SireMol::MolStructureEditor & ( ::SireMol::MolStructureEditor::*removeAllSegments_function_type)(  ) ;
            removeAllSegments_function_type removeAllSegments_function_value( &::SireMol::MolStructureEditor::removeAllSegments );
            
            MolStructureEditor_exposer.def( 
                "removeAllSegments"
                , removeAllSegments_function_value
                , bp::return_self< >()
                , "Remove all segments from this molecule" );
        
        }
        { //::SireMol::MolStructureEditor::rename
        
            typedef ::SireMol::MolStructureEditor & ( ::SireMol::MolStructureEditor::*rename_function_type)( ::SireMol::MolName const & ) ;
            rename_function_type rename_function_value( &::SireMol::MolStructureEditor::rename );
            
            MolStructureEditor_exposer.def( 
                "rename"
                , rename_function_value
                , ( bp::arg("newname") )
                , bp::return_self< >()
                , "Rename this molecule to newname" );
        
        }
        { //::SireMol::MolStructureEditor::renumber
        
            typedef ::SireMol::MolStructureEditor & ( ::SireMol::MolStructureEditor::*renumber_function_type)( ::SireMol::MolNum ) ;
            renumber_function_type renumber_function_value( &::SireMol::MolStructureEditor::renumber );
            
            MolStructureEditor_exposer.def( 
                "renumber"
                , renumber_function_value
                , ( bp::arg("newnum") )
                , bp::return_self< >()
                , "Renumber this molecule to have the number newnum" );
        
        }
        { //::SireMol::MolStructureEditor::renumber
        
            typedef ::SireMol::MolStructureEditor & ( ::SireMol::MolStructureEditor::*renumber_function_type)(  ) ;
            renumber_function_type renumber_function_value( &::SireMol::MolStructureEditor::renumber );
            
            MolStructureEditor_exposer.def( 
                "renumber"
                , renumber_function_value
                , bp::return_self< >()
                , "Give this molecule a new, unique ID number" );
        
        }
        { //::SireMol::MolStructureEditor::residue
        
            typedef ::SireMol::ResStructureEditor ( ::SireMol::MolStructureEditor::*residue_function_type)( ::SireMol::ResID const & ) ;
            residue_function_type residue_function_value( &::SireMol::MolStructureEditor::residue );
            
            MolStructureEditor_exposer.def( 
                "residue"
                , residue_function_value
                , ( bp::arg("resid") )
                , bp::release_gil_policy()
                , "Return an editor for the residue at ID resid\nThrow: SireMol::missing_residue\nThrow: SireMol::duplicate_residue\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMol::MolStructureEditor::segment
        
            typedef ::SireMol::SegStructureEditor ( ::SireMol::MolStructureEditor::*segment_function_type)( ::SireMol::SegID const & ) ;
            segment_function_type segment_function_value( &::SireMol::MolStructureEditor::segment );
            
            MolStructureEditor_exposer.def( 
                "segment"
                , segment_function_value
                , ( bp::arg("segid") )
                , bp::release_gil_policy()
                , "Return an editor for the segment at ID segid\nThrow: SireMol::missing_segment\nThrow: SireMol::duplicate_segment\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMol::MolStructureEditor::select
        
            typedef ::SireMol::AtomStructureEditor ( ::SireMol::MolStructureEditor::*select_function_type)( ::SireMol::AtomID const & ) ;
            select_function_type select_function_value( &::SireMol::MolStructureEditor::select );
            
            MolStructureEditor_exposer.def( 
                "select"
                , select_function_value
                , ( bp::arg("atomid") )
                , bp::release_gil_policy()
                , "Return an editor for the atom at ID atomid\nThrow: SireMol::missing_atom\nThrow: SireMol::duplicate_atom\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMol::MolStructureEditor::select
        
            typedef ::SireMol::CGStructureEditor ( ::SireMol::MolStructureEditor::*select_function_type)( ::SireMol::CGID const & ) ;
            select_function_type select_function_value( &::SireMol::MolStructureEditor::select );
            
            MolStructureEditor_exposer.def( 
                "select"
                , select_function_value
                , ( bp::arg("cgid") )
                , bp::release_gil_policy()
                , "Return an editor for the CutGroup at ID cgid\nThrow: SireMol::missing_cutgroup\nThrow: SireMol::duplicate_cutgroup\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMol::MolStructureEditor::select
        
            typedef ::SireMol::ResStructureEditor ( ::SireMol::MolStructureEditor::*select_function_type)( ::SireMol::ResID const & ) ;
            select_function_type select_function_value( &::SireMol::MolStructureEditor::select );
            
            MolStructureEditor_exposer.def( 
                "select"
                , select_function_value
                , ( bp::arg("resid") )
                , bp::release_gil_policy()
                , "Return an editor for the residue at ID resid\nThrow: SireMol::missing_residue\nThrow: SireMol::duplicate_residue\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMol::MolStructureEditor::select
        
            typedef ::SireMol::ChainStructureEditor ( ::SireMol::MolStructureEditor::*select_function_type)( ::SireMol::ChainID const & ) ;
            select_function_type select_function_value( &::SireMol::MolStructureEditor::select );
            
            MolStructureEditor_exposer.def( 
                "select"
                , select_function_value
                , ( bp::arg("chainid") )
                , bp::release_gil_policy()
                , "Return an editor for the chain at ID chainid\nThrow: SireMol::missing_chain\nThrow: SireMol::duplicate_chain\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMol::MolStructureEditor::select
        
            typedef ::SireMol::SegStructureEditor ( ::SireMol::MolStructureEditor::*select_function_type)( ::SireMol::SegID const & ) ;
            select_function_type select_function_value( &::SireMol::MolStructureEditor::select );
            
            MolStructureEditor_exposer.def( 
                "select"
                , select_function_value
                , ( bp::arg("segid") )
                , bp::release_gil_policy()
                , "Return an editor for the segment at ID segid\nThrow: SireMol::missing_segment\nThrow: SireMol::duplicate_segment\nThrow: SireError::invalid_index\n" );
        
        }
        { //::SireMol::MolStructureEditor::selectedAll
        
            typedef bool ( ::SireMol::MolStructureEditor::*selectedAll_function_type)(  ) const;
            selectedAll_function_type selectedAll_function_value( &::SireMol::MolStructureEditor::selectedAll );
            
            MolStructureEditor_exposer.def( 
                "selectedAll"
                , selectedAll_function_value
                , bp::release_gil_policy()
                , "Return whether or not this is a complete molecule" );
        
        }
        { //::SireMol::MolStructureEditor::toString
        
            typedef ::QString ( ::SireMol::MolStructureEditor::*toString_function_type)(  ) const;
            toString_function_type toString_function_value( &::SireMol::MolStructureEditor::toString );
            
            MolStructureEditor_exposer.def( 
                "toString"
                , toString_function_value
                , bp::release_gil_policy()
                , "Return a string representation of this editor" );
        
        }
        { //::SireMol::MolStructureEditor::typeName
        
            typedef char const * ( *typeName_function_type )(  );
            typeName_function_type typeName_function_value( &::SireMol::MolStructureEditor::typeName );
            
            MolStructureEditor_exposer.def( 
                "typeName"
                , typeName_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        { //::SireMol::MolStructureEditor::what
        
            typedef char const * ( ::SireMol::MolStructureEditor::*what_function_type)(  ) const;
            what_function_type what_function_value( &::SireMol::MolStructureEditor::what );
            
            MolStructureEditor_exposer.def( 
                "what"
                , what_function_value
                , bp::release_gil_policy()
                , "" );
        
        }
        MolStructureEditor_exposer.staticmethod( "typeName" );
        MolStructureEditor_exposer.def( "__copy__", &__copy__);
        MolStructureEditor_exposer.def( "__deepcopy__", &__copy__);
        MolStructureEditor_exposer.def( "clone", &__copy__);
        MolStructureEditor_exposer.def( "__rlshift__", &__rlshift__QDataStream< ::SireMol::MolStructureEditor >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        MolStructureEditor_exposer.def( "__rrshift__", &__rrshift__QDataStream< ::SireMol::MolStructureEditor >,
                            bp::return_internal_reference<1, bp::with_custodian_and_ward<1,2> >() );
        MolStructureEditor_exposer.def_pickle(sire_pickle_suite< ::SireMol::MolStructureEditor >());
        MolStructureEditor_exposer.def( "__str__", &__str__< ::SireMol::MolStructureEditor > );
        MolStructureEditor_exposer.def( "__repr__", &__str__< ::SireMol::MolStructureEditor > );
    }

}
