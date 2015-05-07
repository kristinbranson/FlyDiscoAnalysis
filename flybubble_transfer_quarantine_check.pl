#!/usr/bin/perl

use strict;
use Mail::Sendmail;
use File::Find;
use File::Copy;
use File::Path;
use DBI;
use LWP::UserAgent;
use vars qw($DEBUG);

require '/groups/flyprojects/home/olympiad/bin/flyolympiad_shared_functions.pl';
$DEBUG = 0;

my $dbh = connect_to_sage("/groups/flyprojects/home/olympiad/config/SAGE-prod.config");

my $hr_qc_cvterm_id_lookup = get_qc_cvterm_id($dbh);

foreach my $key (sort keys %$hr_qc_cvterm_id_lookup) {
    print "$key $$hr_qc_cvterm_id_lookup{$key}\n";
}

my $browser = LWP::UserAgent->new;

my $term_id =  get_cv_term_id($dbh,"fly_olympiad","not_applicable");

my $incoming_dir = "/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/00_incoming/";
my $quarantine_dir = "/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/00_quarantine/";
my $bowldata_dir = "/groups/sciserv/flyolympiad/Olympiad_Screen/fly_bowl/bowl_data/";

my %expdirs = ();

my $session_id = "NULL";
my $phase_id = "NULL";

my $skip_jira = 0;

my $message = "";

opendir(INCOMING,"$incoming_dir");
while (my $expcontent = readdir(INCOMING)) {
    $skip_jira = 0;
    my $experror_message = "";
    chomp($expcontent);
    next if ($expcontent =~ /^\./);
    
    my $dir_path = $incoming_dir . $expcontent;
    my $new_path = $quarantine_dir . $expcontent;
    my $bowldata_path = $bowldata_dir . $expcontent;
    
    print "Checking $dir_path\n";

    my $has_ufmf = 0;
    my $has_log = 0;
    my $has_metadata = 0;
    my $has_protocol = 0;
    my $has_ufmf_diag = 0;
    my $has_aborted = 0;

    if ($expcontent =~ /notstarted/) {
    	$skip_jira = 1;
    	$experror_message .= "Experiment $expcontent has been flagged as not started.";
    } else {
    	opendir(EXPDIR,$dir_path);
    	while (my $expfoldercontents = readdir(EXPDIR)) {
    	    #print "$expfoldercontents\n" if ($DEBUG);
    	    if ($expfoldercontents =~ /movie.ufmf$/) {
        		my $filepath = $dir_path . "/" . $expfoldercontents;
        		$has_ufmf++ if (-s $filepath);
    	    }
    	    if ($expfoldercontents =~ /Log.txt/) {
        		my $filepath = $dir_path . "/" . $expfoldercontents;
        		$has_log++ if (-s $filepath);
    	    }
            if ($expfoldercontents =~ /Metadata.xml$/) {
        		my $filepath = $dir_path . "/" . $expfoldercontents;
        		$has_metadata++ if (-s $filepath);
    	    }
            if ($expfoldercontents =~ /protocol.mat/) {
        		my $filepath = $dir_path . "/" . $expfoldercontents;
        		$has_protocol++ if (-s $filepath);
    	    }
            if ($expfoldercontents =~ /ufmf_diagnostics.txt/) {
        		my $filepath = $dir_path . "/" . $expfoldercontents;
        		$has_ufmf_diag++ if (-s $filepath);
    	    }
    	    if ($expfoldercontents =~ /ABORTED$/) {
    		    #$skip_jira = 1;
                my $filepath = $dir_path . "/" . $expfoldercontents;
                $has_aborted++ if (-s $filepath);
            }
    	}
    	closedir(EXPDIR);
    }
    
    if ($DEBUG) { print "$has_ufmf $has_log $has_metadata $has_protocol $has_ufmf_diag\n"; }
    if ($has_aborted) {
        $experror_message .= "Fly bowl experiment has been aborted\n";
    }
    unless ($has_ufmf) {
	    $experror_message .= "movie.ufmf is missing or empty\n";
    }
    unless ($has_log) {
        $experror_message .= "Log.txt is missing or empty\n";
    }
    unless ($has_metadata) {
        $experror_message .= "Metadata.xml is missing or empty\n";
    }
    unless ($has_temperature) {
        $experror_message .= "protocol.mat is missing or empty\n";
    }
    unless ($has_ufmf_diag) {
        $experror_message .= "ufmf_diagnostics.txt is missing or empty\n";
    }


    my $sage_exp_id;

    if ($has_metadata) {
    	# check if experiment already loaded
    	$sage_exp_id = get_experiment_id($dbh,$expcontent);
    	if ($sage_exp_id) {
    	    print "$expcontent already in sage\n";
    	} else {
    	    #Load MetaData;
    	    print "loading $expcontent\n";
    	    my $cmd = "/groups/flyprojects/home/olympiad/bowl_bin/wrapper_FlyBowlStoreExperimentProd2.sh $dir_path";
    	    print "$cmd\n";
    	    unless ($DEBUG) {
        		system($cmd);
        		$sage_exp_id = get_experiment_id($dbh,$expcontent);
        		my $exp_name = get_experiment_name($dbh, $sage_exp_id);

        		#my $flybowl_exp_name = "FlyBowl_" . $exp_name;
        		#my $update_name = "update experiment set name = '$flybowl_exp_name' where id = $sage_exp_id";
        		#run_mod($dbh,$update_name);
    	    }
    	}
    }

    my $automated_pf_category = "";

    if ($sage_exp_id) {
	    print "EXP_ID: $sage_exp_id\n" if ($DEBUG);
        my $screen_type = get_screen_type($dbh,$sage_exp_id);
    	my $analysis_type = "";
    	if ($screen_type =~ /non_olympiad/) {
    	    $analysis_type = "current_" . $screen_type;
    	}
        # check to see if we should fail experiment by analyzing metadata
        system("/groups/flyprojects/home/olympiad/bowl_bin/wrapper_FlyBowlAutomaticChecks_Incoming.sh $dir_path $analysis_type");
    	#make sure every loaded experiments gets an automated_pf flag.
    	insert_exp_automated_pf($dbh,$sage_exp_id,"U");
    	# load automated checks value for automated_pf
    	unless (-e "$dir_path/automatic_checks_incoming_results.txt") {
    	    $experror_message .= "automated_checks_incoming_results failed to generate result.\n";
    	}
    	system("/groups/flyprojects/home/olympiad/bowl_bin/flybowl_load_automated_checks.pl $dir_path automatic_checks_incoming_results.txt");
    	system("/groups/flyprojects/home/olympiad/bowl_bin/flybowl_inline_dataloader.pl $dir_path QuickStats.txt");
    	system("/groups/flyprojects/home/olympiad/bowl_bin/flybowl_ufmf_diagnostics_dataloader.pl $dir_path ufmf_diagnostics.txt");
    	system("/groups/flyprojects/home/olympiad/bowl_bin/flybowl_temperature_stream.pl $dir_path");
    	
    	my $fsp_cv = "fly_olympiad";
    	my $fsp_cvterm = "file_system_path";
    	load_experiment_property($dbh,$sage_exp_id,$fsp_cv,$fsp_cvterm,$bowldata_path);
    } else {
    	$automated_pf_category = "sage_load_failed";
    	$experror_message .= "Metadata.xml failed to load into SAGE.\n";
    }

    # check automated pf in database if F error experiment, if P then ok it for analysis
    my ($prop_id, $automated_pf) = get_experiment_property_val($dbh, $sage_exp_id, "fly_olympiad_qc", "automated_pf");
    my $apfc_id = 0;
    if ($automated_pf eq "F") {
        # Here
	    ($apfc_id, $automated_pf_category) = get_experiment_property_val($dbh, $sage_exp_id, "fly_olympiad_qc", "automated_pf_category");
        my ($nc_id, $notes_cure) = get_experiment_property_val($dbh, $sage_exp_id, "fly_olympiad_qc", "notes_curation");
	    $experror_message .= "Error: automated checks for incoming experiments flagged experiment for Fail\n$notes_cure\n";
    }

    if ($DEBUG) { print "$experror_message\n"; }
    
    if ($experror_message) {
    	print "Detected Errors:\n$experror_message";
    	$experror_message .= "File System Path: $new_path\n";
    	unless($DEBUG) {
    	    move ("$dir_path", "$new_path" );
    	}

    	# set automated_pf to fail
    	if ($sage_exp_id) {
    	    insert_exp_automated_pf($dbh,$sage_exp_id,"F");
    	    insert_exp_manual_pf($dbh,$sage_exp_id,"U");
    	}

    	my $label = "uncategorized_transfer_error";

    	if ($automated_pf_category) {
    	    $label = "$automated_pf_category";
    	}

    	#jira ticket submission	
    	my %jira_ticket_params;
    	$jira_ticket_params{'lwp_handle'} = $browser;
    	$jira_ticket_params{'jira_project_pid'} = 10102;
    	$jira_ticket_params{'issue_type_id'} = 6;
    	$jira_ticket_params{'labels'} = $label;
    	$jira_ticket_params{'summary'} = "Fly Bowl Transfer Error Detected $expcontent";
    	$jira_ticket_params{'description'} = $experror_message;
    	$jira_ticket_params{'file_path'} = $new_path;
    	$jira_ticket_params{'error_type'} = "";
    	unless($skip_jira) {
    	    print "Errors found submitting Jira Ticket\n";
    	    submit_jira_ticket(\%jira_ticket_params);
    	}
    	$message .= $experror_message;
    } else {
    	if ($sage_exp_id) {
    	    #insert_exp_automated_pf($dbh,$sage_exp_id,"P");
    	    insert_exp_manual_pf($dbh,$sage_exp_id,"U");
    	}
    }

}

#mail notification
closedir(INCOMING);

$dbh->disconnect();


if ($message) {
    $message = "Fly Bowl data transfer check\n" . $message;
    my $subject = "[Olympiad Fly Bowl Data Transfer]Quarantine experiments that are incomplete";
    
    #send_email('@janelia.hhmi.org','olympiad@janelia.hhmi.org', $subject, $message,'midgleyf@janelia.hhmi.org');
    send_email('midgleyf@janelia.hhmi.org','olympiad@janelia.hhmi.org', $subject, $message);
}

exit;



sub get_tube_avi_count {
    my ($temperature_dir) = @_;
    my $dir_avi_count = 0;
    opendir(TEMPERDIR,"$temperature_dir");
    while (my $avi = readdir(TEMPERDIR)) {
        if ($avi =~ /seq\d+_tube\d+\.avi/) {
            my $avi_path = $temperature_dir . "/" . $avi;
            my $file_size = 0;
            $file_size = -s "$avi_path";
            #print "$avi_path $file_size\n";
            if ($file_size > 0) {
        	$dir_avi_count++;
            }
        }
    }
    return($dir_avi_count);
}

sub get_screen_type {
    my($dbh,$sage_exp_id) = @_;
    my $sql = qq~
select p.value 
from experiment e, experiment_property p, cv_term c 
where e.id = $sage_exp_id  
and e.id = p.experiment_id 
and p.type_id = c.id 
and c.name = 'screen_type' 
~;
    my @toarray = do_sql($dbh, $sql);
    return($toarray[0]);
}

sub get_qc_cvterm_id {
    my($dbh) = @_;
    my %qc_cvterm_ids;

    my $qc_cv_sql = "select ct.name, ct.id from cv_term ct, cv c where c.name = 'fly_olympiad_qc' and ct.cv_id = c.id";

    my @qc_rows = do_sql($dbh,$qc_cv_sql);

    foreach my $row (@qc_rows) {
        my ($termname,$termid) = split(/\t/,$row);

        #print "$termname,$termid\n";
        $qc_cvterm_ids{$termname} = $termid;
    }
    return(\%qc_cvterm_ids);
}
