--schema for regulatory information collected from other sources

create schema harbison2004;

set search_path to harbison2004;

create table readme(
	id	serial primary key,
	name	varchar,
	description	varchar
	);

--matrix is actually factor, for backward compatibility
create table matrix(
	id	serial primary key,
	mt_acc	varchar unique,
	mt_id	varchar unique,
	tax_id	integer,	--12-09-05
	tf_name	varchar,
	tf_description	varchar,
	factors	varchar[],
	pwm	float[][],
	consensus	varchar,
	basis	varchar,
	sites	varchar[],
	site_in_matrix_accs	varchar[],
	comment	varchar,
	reference	varchar[],
	entrezgene_id	integer	--12-09-05 new id
	);
	
create table binding_site(
	id	serial primary key,
	mt_id	varchar references matrix(mt_id) on delete cascade on update cascade,
	chromosome	varchar,
	strand	varchar(1),
	bs_genome_start	integer,
	bs_genome_end	integer,
	bs_disp_start	integer,
	bs_disp_end	integer,
	prom_id	integer not null,	--the gene id
	tax_id	integer,	--12-09-05
	core_similarity_score	float,
	matrix_similarity_score	float,
	sequence	varchar,
	comment	varchar	--condition, etc.
	);


CREATE OR REPLACE FUNCTION before_binding_site() RETURNS trigger AS '
	DECLARE
		prom_seq_position RECORD;
	BEGIN
		IF NEW.mt_id='''' THEN
			NEW.mt_id := null;
		END IF;
		IF NEW.chromosome='''' THEN
			NEW.chromosome := null;
		END IF;
		IF NEW.strand='''' THEN
			NEW.strand := null;
		END IF;
		IF NEW.bs_genome_start < 0 THEN
			NEW.bs_genome_start := null;
		END IF;
		IF NEW.bs_genome_end < 0 THEN
			NEW.bs_genome_end := null;
		END IF;
		IF NEW.tax_id < 0 THEN
			NEW.tax_id := null;
		END IF;
		IF NEW.core_similarity_score < 0 THEN
			NEW.core_similarity_score := null;
		END IF;
		IF NEW.matrix_similarity_score < 0 THEN
			NEW.matrix_similarity_score := null;
		END IF;
		IF NEW.sequence = '''' THEN
			NEW.sequence := null;
		END IF;
		IF NEW.comment = '''' THEN
			NEW.comment := null;
		END IF;
		RETURN NEW;
	END;
	' LANGUAGE plpgsql;

CREATE TRIGGER before_binding_site BEFORE INSERT OR UPDATE ON binding_site
	FOR EACH ROW EXECUTE PROCEDURE before_binding_site();
