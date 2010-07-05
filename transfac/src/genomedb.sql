-- 2008-07-27 this is deprecated. GenomeDB.py defines tables and etc. genomedb_tablenames.txt contains names of all tables. pymodule/OutputSQLTrigger.py reads genomedb_tablenames.txt to output triggers in sql syntax.
-- 2008-07-06 add foreign keys and schema name change: sequence -> genome --
-- 11-12-05 genomedb --

create schema genome;
set search_path to genome;
create table sequence_type(
	id integer primary key,
	type varchar
	);

create table raw_sequence(
	id serial primary key,
	acc_ver varchar,
	start integer,
	stop integer,
	sequence varchar(10000)	--the unit is 1kb--
	);


create table annot_assembly(
	id serial,
	gi integer primary key,
	acc_ver varchar unique,
	acc varchar,
	version integer,
	tax_id integer,
	chromosome varchar,
	start integer,
	stop integer,
	orientation varchar(1),
	sequence varchar,
	raw_sequence_start_id integer references raw_sequence (id) on delete restrict on update cascade,
	seq_type integer references sequence_type (id) on delete restrict on update cascade,
	comment varchar
	);

create table entrezgene_mapping(
	id serial,
	gene_id varchar primary key,
	tax_id integer,
	genomic_acc_ver varchar references annot_assembly (acc_ver) on delete restrict on update cascade,
	genomic_gi integer,
	strand varchar,
	start integer,
	stop integer,
	mrna_acc_ver varchar,
	mrna_gi integer,
	mrna_start integer[],
	mrna_stop integer[],
	cds_acc_ver varchar,
	cds_gi integer,
	cds_start integer[],
	cds_stop integer[],
	comment varchar
	);

create table gene(
	id serial,
	tax_id integer,
	gene_id integer primary key,
	gene_symbol varchar,
	locustag varchar,
	synonyms varchar,
	dbxrefs varchar,
	chromosome varchar,
	map_location varchar,
	description varchar,
	type_of_gene varchar,
	symbol_from_nomenclature_authority varchar,
	full_name_from_nomenclature_authority varchar,
	nomenclature_status varchar,
	other_designations varchar,
	modification_date date
);

CREATE TABLE gene2go(
	id serial primary key,
    tax_id integer,
    gene_id integer references gene (gene_id) on delete restrict on update cascade,
    go_id character varying,
    evidence character varying,
    go_qualifier character varying,
    go_description character varying,
    pubmed_ids character varying,
    category varchar
);

CREATE TABLE gene_symbol2id (
	id serial primary key,
    tax_id integer,
    gene_symbol character varying,
    gene_id integer references gene (gene_id) on delete restrict on update cascade,
    symbol_type character varying
);

create table readme(
        id      serial primary key,
        title   varchar,
        description     varchar,
        created_by        varchar default current_user,
        updated_by      varchar,
        date_created    timestamp WITH TIME ZONE default current_timestamp,
        date_updated   TIMESTAMP WITH TIME ZONE
        );

create function insert_readme() returns trigger as '
        begin
        new.created_by := current_user;
        return new;
        end;
        '
        language plpgsql;

create trigger insert_readme before insert on readme
        for each row execute procedure insert_readme();

create function update_readme() returns trigger as '
        begin
        new.date_updated := current_timestamp;
        new.updated_by := current_user;
        return new;
        end;
        '
        language plpgsql;

create trigger update_readme before update on readme
        for each row execute procedure update_readme();

insert into sequence_type(type) values('chromosome sequence');
insert into sequence_type(type) values('CDS');
insert into sequence_type(type) values('exon');
insert into sequence_type(type) values('intron');
insert into sequence_type(type) values("5' UTR");
insert into sequence_type(type) values("3' UTR");
insert into sequence_type(type) values("other");