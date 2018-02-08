# Copyright 2011 Elhanan Borenstein and Rogan Carr
#  Contact: rogan@uw.edu
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

package NetSeed;
require Exporter;
use strict;

our @ISA = ("Exporter");
our @EXPORT = qw( CalculateSeeds Seeds GroupedSeeds NonSeeds IgnoredNodes AllNodes AllEdges NumElements );
our @EXPORT_OK = qw();

# Globals (but not exported)
#  These are accessed externally through function calls.
our %AllNodes =();
our %AllEdges =();
our %v_nodes_pruned;
our %v_seeds = ();
our %v_seeds_grouped = ();
our %v_non_seeds = ();
our $totalCount;

sub CalculateSeeds {
  
  my $netRef = shift;
  my $SEPARATOR = shift;
  my $GIANT_COMP = shift;
  my $MINIMAL_COMP_SIZE = shift;
  
  # HASHES INIT
  # ===========  
  my %v_nodes = ();
  my %v_nodes_order = ();
  my $nof_nodes = 0;
  %v_nodes_pruned = ();
  my %gEdges;
  my %gRevEdges;

  $totalCount = 0;

  # ============================================
  # ==== READ NETWORK
  # ============================================

  #   open (RL_E_INFILE, $NETWORK_FILE) or die "Kozyol  pp !\n";
  my $edgeCount = 0;
  while (<$netRef>) {
    my $line_e=$_; chomp($line_e);
    my $edge = $line_e;
    
    if (length($edge) >= 3) # A VALID LINE CANNOT BE SHORTER THAN THIS
      {
	my ($edge_in, $edge_out) = split(/$SEPARATOR/, $edge);
	
	$gEdges{$edge_in}->{$edge_out} = 1;
	$gRevEdges{$edge_out}->{$edge_in} = 1;

	# We need a starting copy to draw later
	$AllEdges{$edge_in}->{$edge_out} = 1;

	# KEEP TRACK OF NODES
	if (not $v_nodes{$edge_in}) {
	  $nof_nodes++;  $v_nodes{$edge_in} = 1;  $v_nodes_order{$nof_nodes} = $edge_in;
	}

	if (not $v_nodes{$edge_out}) {
	  $nof_nodes++;  $v_nodes{$edge_out} = 1;  $v_nodes_order{$nof_nodes} = $edge_out;
	}
	$edgeCount++;
      }
  }

  %AllNodes = %v_nodes;

  $totalCount = $edgeCount + $nof_nodes;

  # ============================================
  # ==== TRIM SMALL COMPONENTS
  # ============================================
    
  # CALL DFS
  my ($node_finish_order, $component, $nc) = DFS(\%gEdges,\%v_nodes_order,$nof_nodes,\%gRevEdges);
  my %node_component = %{$nc};
	
  # CALCULATE THE SIZE OF EACH CONNECTED COMPONENT AND FIND THE LARGEST COMPONENT
  my %v_comp_size = ();
  my $largest_component = -1;
  my $largest_component_size = 0;
  for (my $index=1; $index<=$nof_nodes; $index++) {
    my $node = $v_nodes_order{$index};
    my $curr_component = $node_component{$node};
    if ($v_comp_size{$curr_component}) {
      $v_comp_size{$curr_component}++;
    } else {
      $v_comp_size{$curr_component}=1;
    }
		
    if ($v_comp_size{$curr_component} > $largest_component_size) {
      $largest_component_size = $v_comp_size{$curr_component};
      $largest_component = $curr_component;
    }
  }

  # DELETE EDGES OF TOO SMALL COMPONENTS OR ANYTHING BUT THE LARGEST COMPONENT
  my $comp_TH;
  if ($GIANT_COMP == 1) {
    $comp_TH = $largest_component_size;
  } else {
    $comp_TH = $MINIMAL_COMP_SIZE;
  }

  foreach my $edge_in (keys %gEdges) {
    foreach my $edge_out (keys %{$gEdges{$edge_in}}) {
      if ($v_comp_size{$node_component{$edge_in}} < $comp_TH) {
	delete $gEdges{$edge_in}->{$edge_out};
	delete $gRevEdges{$edge_out}->{$edge_in};
	
	$v_nodes_pruned{$edge_in} = 1;
	$v_nodes_pruned{$edge_out} = 1;
      }
    }
  }

  # COUNT NODES AGAIN AND MAKE LISTS (AFTER TRIMMING)
  %v_nodes = ();
  %v_nodes_order = ();
  $nof_nodes = 0;

  foreach my $edge_in (keys %gEdges) {
    foreach my $edge_out (keys %{$gEdges{$edge_in}}) {

      # KEEP TRACK OF NODES
      if (not $v_nodes{$edge_in}) {
	$nof_nodes++;  $v_nodes{$edge_in} = 1;  $v_nodes_order{$nof_nodes} = $edge_in;
      }
      
      if (not $v_nodes{$edge_out}) {
	$nof_nodes++;  $v_nodes{$edge_out} = 1;  $v_nodes_order{$nof_nodes} = $edge_out;
      }
    }
  }
  
  # Network created; Find the seeds

  my $MINIMAL_SEED_PROB = 0; 
  
  %v_seeds = ();
  %v_seeds_grouped = ();
  %v_non_seeds = ();
  
  # ===========================================
  # ==== RUN DFS TWICE TO FIND SCC
  # ===========================================

  my $node_order_dfs;
  ($node_order_dfs,$component,$nc) = DFS(\%gEdges, \%v_nodes_order,$nof_nodes);
  
  ($node_order_dfs,$component,$nc) = DFS(\%gRevEdges, \%{$node_order_dfs},$nof_nodes);
  %node_component = %{$nc};
  
  my $nof_components = $component;
  
  # ===============================================
  # ==== CREATE SCC COMPONENT AND IDENTIFY SOURCES
  # ===============================================
  
  my %v_scc_group = ();
  my %v_scc_group_size = ();
  my %v_scc_group_source = ();
  
  for (my $index=1; $index<=$nof_nodes; $index++) {
    my $node = $v_nodes_order{$index};
    my $curr_component = $node_component{$node};
    
    if ($v_scc_group{$curr_component}) {
      $v_scc_group{$curr_component} .= "$SEPARATOR$node"; $v_scc_group_size{$curr_component}++;
    } else {
      $v_scc_group{$curr_component}  = "$node"; $v_scc_group_size{$curr_component} = 1; $v_scc_group_source{$curr_component} = 1;
    } 
    # Initially all SCC are considered source
  }
  
  # MARK ALL SCC THAT ARE NOT SOURCE SCCs
  foreach my $edge_in (keys %gEdges) {
    foreach my $edge_out (keys %{$gEdges{$edge_in}}) {
      # Chech if edge is between SCC and mark the SCC in which edge_out is located as NOT a source SCC
      if ($node_component{$edge_out} ne $node_component{$edge_in}) {
	$v_scc_group_source{$node_component{$edge_out}} = 0;
      } 
    }
  }
  
  # ===========================================
  # ==== STORE SEED INFO
  # ===========================================
  
  foreach my $node (sort(keys %node_component)) { 
    my $node_seed_prob = $v_scc_group_source{$node_component{$node}} / $v_scc_group_size{$node_component{$node}}; 
    
    if ($node_seed_prob > $MINIMAL_SEED_PROB) {
      $v_seeds{$node} = $node_seed_prob;
    } else {
      $v_non_seeds{$node} = 1;
    }
  }
  
  for (my $index=1; $index<=$nof_components; $index++) {
    if ($v_scc_group_source{$index} == 1) {
      $v_seeds_grouped{$v_scc_group{$index}} = 1;
    }
  }

}

sub Seeds {
  return %v_seeds;
}

sub GroupedSeeds {
  return %v_seeds_grouped;
}

sub NonSeeds {
  return %v_non_seeds;
}

sub IgnoredNodes {
  return %v_nodes_pruned;
}

sub AllNodes {
  return %AllNodes;
}

sub AllEdges {
  return %AllEdges;
}

sub NumElements {
  return $totalCount;
}

# ======================================================================
# ======================================================================
# ======================================================================

sub DFS {
  my $Edges = shift;
  my $nodeOrder = shift;
  my $nof_nodes = shift;
  my $revEdges = shift;
    
  my @L = ();
  my $finish = 0;
  my %node_visit = ();
  my %node_finish = ();
  my %node_finish_order = ();
  
  my $component = 0;
  my %node_component = ();
  
  for (my $index=$nof_nodes; $index>=1; $index--) {
    my $x = $nodeOrder->{$index};
    if ($node_visit{$x}) {
      next;
    }
    
    $component++;
    
    $node_component{$x} = $component;
    push(@L,$x);
    $node_visit{$x} = 1;
    
    
    SEARCH($x,\%node_visit,$Edges,\@L);
    if (defined $revEdges) {
      SEARCH($x,\%node_visit,$revEdges,\@L);
    }
    
    while ($#L+1 > 0) {
      my $w = pop(@L);
      if ($node_visit{$w}) {
	if (not $node_finish{$w}) {
	  $node_finish{$w} = 1; $finish++; $node_finish_order{$finish} = $w;
	}
      } else {
	$node_component{$w} = $component;
	push(@L,$w);
	$node_visit{$w} = 1;
	
	SEARCH($w,\%node_visit,$Edges,\@L);
	if (defined $revEdges) {
	  SEARCH($w,\%node_visit,$revEdges,\@L);
	}
	
      }
    }
  }
  return (\%node_finish_order, $component, \%node_component)
}

sub SEARCH {
  my $v = shift;
  my $visitedNodes = shift;
  my $Edges = shift;
  my $L = shift;

  if (defined $Edges->{$v}) {
    foreach my $outEdge (keys %{$Edges->{$v}}) {
      unless ( $visitedNodes->{$outEdge} ) {
	push(@{$L}, $outEdge);
      }
    }
  }  
}


1;
