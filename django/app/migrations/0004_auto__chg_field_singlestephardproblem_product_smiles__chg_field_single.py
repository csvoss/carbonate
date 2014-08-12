# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):

        # Changing field 'SingleStepHardProblem.product_smiles'
        db.alter_column(u'app_singlestephardproblem', 'product_smiles', self.gf('django.db.models.fields.TextField')())

        # Changing field 'SingleStepHardProblem.reactant_smiles'
        db.alter_column(u'app_singlestephardproblem', 'reactant_smiles', self.gf('django.db.models.fields.TextField')())

        # Changing field 'SingleStepProblem.product_smiles'
        db.alter_column(u'app_singlestepproblem', 'product_smiles', self.gf('django.db.models.fields.TextField')())

        # Changing field 'SingleStepProblem.reactant_smiles'
        db.alter_column(u'app_singlestepproblem', 'reactant_smiles', self.gf('django.db.models.fields.TextField')())

        # Changing field 'PredictProductsProblem.reactant_smiles'
        db.alter_column(u'app_predictproductsproblem', 'reactant_smiles', self.gf('django.db.models.fields.TextField')())

        # Changing field 'PredictProductsProblem.correct_answer'
        db.alter_column(u'app_predictproductsproblem', 'correct_answer', self.gf('django.db.models.fields.TextField')())

        # Changing field 'PredictProductsProblem.incorrect_answer3'
        db.alter_column(u'app_predictproductsproblem', 'incorrect_answer3', self.gf('django.db.models.fields.TextField')())

        # Changing field 'PredictProductsProblem.incorrect_answer2'
        db.alter_column(u'app_predictproductsproblem', 'incorrect_answer2', self.gf('django.db.models.fields.TextField')())

        # Changing field 'PredictProductsProblem.incorrect_answer1'
        db.alter_column(u'app_predictproductsproblem', 'incorrect_answer1', self.gf('django.db.models.fields.TextField')())

        # Changing field 'Synthesis.product_smiles'
        db.alter_column(u'app_synthesis', 'product_smiles', self.gf('django.db.models.fields.TextField')())

        # Changing field 'Synthesis.reactant_smiles'
        db.alter_column(u'app_synthesis', 'reactant_smiles', self.gf('django.db.models.fields.TextField')())

    def backwards(self, orm):

        # Changing field 'SingleStepHardProblem.product_smiles'
        db.alter_column(u'app_singlestephardproblem', 'product_smiles', self.gf('django.db.models.fields.CharField')(max_length=100))

        # Changing field 'SingleStepHardProblem.reactant_smiles'
        db.alter_column(u'app_singlestephardproblem', 'reactant_smiles', self.gf('django.db.models.fields.CharField')(max_length=100))

        # Changing field 'SingleStepProblem.product_smiles'
        db.alter_column(u'app_singlestepproblem', 'product_smiles', self.gf('django.db.models.fields.CharField')(max_length=100))

        # Changing field 'SingleStepProblem.reactant_smiles'
        db.alter_column(u'app_singlestepproblem', 'reactant_smiles', self.gf('django.db.models.fields.CharField')(max_length=100))

        # Changing field 'PredictProductsProblem.reactant_smiles'
        db.alter_column(u'app_predictproductsproblem', 'reactant_smiles', self.gf('django.db.models.fields.CharField')(max_length=100))

        # Changing field 'PredictProductsProblem.correct_answer'
        db.alter_column(u'app_predictproductsproblem', 'correct_answer', self.gf('django.db.models.fields.CharField')(max_length=100))

        # Changing field 'PredictProductsProblem.incorrect_answer3'
        db.alter_column(u'app_predictproductsproblem', 'incorrect_answer3', self.gf('django.db.models.fields.CharField')(max_length=100))

        # Changing field 'PredictProductsProblem.incorrect_answer2'
        db.alter_column(u'app_predictproductsproblem', 'incorrect_answer2', self.gf('django.db.models.fields.CharField')(max_length=100))

        # Changing field 'PredictProductsProblem.incorrect_answer1'
        db.alter_column(u'app_predictproductsproblem', 'incorrect_answer1', self.gf('django.db.models.fields.CharField')(max_length=100))

        # Changing field 'Synthesis.product_smiles'
        db.alter_column(u'app_synthesis', 'product_smiles', self.gf('django.db.models.fields.CharField')(max_length=100))

        # Changing field 'Synthesis.reactant_smiles'
        db.alter_column(u'app_synthesis', 'reactant_smiles', self.gf('django.db.models.fields.CharField')(max_length=100))

    models = {
        u'api.property': {
            'Meta': {'object_name': 'Property'},
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50'})
        },
        u'api.reagent': {
            'Meta': {'object_name': 'Reagent'},
            'diagram_name': ('django.db.models.fields.CharField', [], {'max_length': '50', 'blank': 'True'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'is_solvent': ('django.db.models.fields.BooleanField', [], {'default': 'False'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'properties': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['api.Property']", 'null': 'True', 'blank': 'True'}),
            'smiles': ('django.db.models.fields.CharField', [], {'max_length': '100', 'blank': 'True'})
        },
        u'api.reagentset': {
            'Meta': {'object_name': 'ReagentSet'},
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'name': ('django.db.models.fields.CharField', [], {'max_length': '50', 'blank': 'True'}),
            'reagents': ('django.db.models.fields.related.ManyToManyField', [], {'related_name': "'reagentSets'", 'symmetrical': 'False', 'to': u"orm['api.Reagent']"}),
            'solvent': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['api.Reagent']", 'null': 'True', 'blank': 'True'}),
            'solvent_properties': ('django.db.models.fields.related.ManyToManyField', [], {'symmetrical': 'False', 'to': u"orm['api.Property']", 'null': 'True', 'blank': 'True'})
        },
        u'app.predictproductsproblem': {
            'Meta': {'object_name': 'PredictProductsProblem'},
            'correct_answer': ('django.db.models.fields.TextField', [], {}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'incorrect_answer1': ('django.db.models.fields.TextField', [], {}),
            'incorrect_answer2': ('django.db.models.fields.TextField', [], {}),
            'incorrect_answer3': ('django.db.models.fields.TextField', [], {}),
            'reactant_smiles': ('django.db.models.fields.TextField', [], {}),
            'reagents': ('django.db.models.fields.CharField', [], {'max_length': '100'})
        },
        u'app.singlestephardproblem': {
            'Meta': {'object_name': 'SingleStepHardProblem'},
            'answer': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['api.ReagentSet']"}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'product_smiles': ('django.db.models.fields.TextField', [], {}),
            'reactant_smiles': ('django.db.models.fields.TextField', [], {})
        },
        u'app.singlestepproblem': {
            'Meta': {'object_name': 'SingleStepProblem'},
            'correct_answer': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'product_smiles': ('django.db.models.fields.TextField', [], {}),
            'reactant_smiles': ('django.db.models.fields.TextField', [], {})
        },
        u'app.synthesis': {
            'Meta': {'object_name': 'Synthesis'},
            'correct_answer': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'product_smiles': ('django.db.models.fields.TextField', [], {}),
            'reactant_smiles': ('django.db.models.fields.TextField', [], {})
        }
    }

    complete_apps = ['app']