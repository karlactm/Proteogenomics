package BLASTProtein;

sub new {
    my $class = shift;
    my $self = {
        _query_id => shift,
        _subject_id => shift,
        _identity => shift,   
        _positives => shift,
        _bit_score => shift,
        _subject_cov => shift,
        _hsp_cov => shift
    };

    bless $self, $class;
    return $self;
}

sub getQueryId {
    my $self = shift;
    return $self->{_query_id};
}

sub getSubjectId {
    my $self = shift;
    return $self->{_subject_id};
}

sub getIdentity {
    my $self = shift;
    return $self->{_identity};
}

sub getPositives {
    my $self = shift;
    return $self->{_positives};
}

sub getBitScore {
    my $self = shift;
    return $self->{_bit_score};
}

sub getSubjectCov {
    my $self = shift;
    return $self->{_subject_cov};
}

sub getHspCov {
    my $self = shift;
    return $self->{_hsp_cov};
}

1;
