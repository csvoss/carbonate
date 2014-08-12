# -*- coding: utf-8 -*-
from south.utils import datetime_utils as datetime
from south.db import db
from south.v2 import SchemaMigration
from django.db import models


class Migration(SchemaMigration):

    def forwards(self, orm):
        # Adding model 'SingleStepHardProblem'
        db.create_table(u'app_singlestephardproblem', (
            ('id', self.gf('django.db.models.fields.AutoField')(primary_key=True)),
            ('reactant_smiles', self.gf('django.db.models.fields.CharField')(max_length=100)),
            ('product_smiles', self.gf('django.db.models.fields.CharField')(max_length=100)),
            ('answer', self.gf('django.db.models.fields.related.ForeignKey')(to=orm['api.ReagentSet'])),
        ))
        db.send_create_signal(u'app', ['SingleStepHardProblem'])

        # Deleting field 'SingleStepProblem.productSmiles'
        db.delete_column(u'app_singlestepproblem', 'productSmiles')

        # Deleting field 'SingleStepProblem.reactantSmiles'
        db.delete_column(u'app_singlestepproblem', 'reactantSmiles')

        # Deleting field 'SingleStepProblem.correctAnswer'
        db.delete_column(u'app_singlestepproblem', 'correctAnswer')

        # Adding field 'SingleStepProblem.reactant_smiles'
        db.add_column(u'app_singlestepproblem', 'reactant_smiles',
                      self.gf('django.db.models.fields.CharField')(default='CCC=CC', max_length=100),
                      keep_default=False)

        # Adding field 'SingleStepProblem.product_smiles'
        db.add_column(u'app_singlestepproblem', 'product_smiles',
                      self.gf('django.db.models.fields.CharField')(default='CCCCC', max_length=100),
                      keep_default=False)

        # Adding field 'SingleStepProblem.correct_answer'
        db.add_column(u'app_singlestepproblem', 'correct_answer',
                      self.gf('django.db.models.fields.CharField')(default='Magic', max_length=100),
                      keep_default=False)

        # Deleting field 'PredictProductsProblem.reactantSmiles'
        db.delete_column(u'app_predictproductsproblem', 'reactantSmiles')

        # Deleting field 'PredictProductsProblem.correctAnswer'
        db.delete_column(u'app_predictproductsproblem', 'correctAnswer')

        # Adding field 'PredictProductsProblem.reactant_smiles'
        db.add_column(u'app_predictproductsproblem', 'reactant_smiles',
                      self.gf('django.db.models.fields.CharField')(default='CCCCC', max_length=100),
                      keep_default=False)

        # Adding field 'PredictProductsProblem.correct_answer'
        db.add_column(u'app_predictproductsproblem', 'correct_answer',
                      self.gf('django.db.models.fields.CharField')(default='Magic', max_length=100),
                      keep_default=False)

        # Adding field 'PredictProductsProblem.incorrect_answer1'
        db.add_column(u'app_predictproductsproblem', 'incorrect_answer1',
                      self.gf('django.db.models.fields.CharField')(default='More Magic', max_length=100),
                      keep_default=False)

        # Adding field 'PredictProductsProblem.incorrect_answer2'
        db.add_column(u'app_predictproductsproblem', 'incorrect_answer2',
                      self.gf('django.db.models.fields.CharField')(default='CCCCl', max_length=100),
                      keep_default=False)

        # Adding field 'PredictProductsProblem.incorrect_answer3'
        db.add_column(u'app_predictproductsproblem', 'incorrect_answer3',
                      self.gf('django.db.models.fields.CharField')(default='CCCCBr', max_length=100),
                      keep_default=False)

        # Deleting field 'Synthesis.productSmiles'
        db.delete_column(u'app_synthesis', 'productSmiles')

        # Deleting field 'Synthesis.reactantSmiles'
        db.delete_column(u'app_synthesis', 'reactantSmiles')

        # Deleting field 'Synthesis.correctAnswer'
        db.delete_column(u'app_synthesis', 'correctAnswer')

        # Adding field 'Synthesis.reactant_smiles'
        db.add_column(u'app_synthesis', 'reactant_smiles',
                      self.gf('django.db.models.fields.CharField')(default='CCCC=C', max_length=100),
                      keep_default=False)

        # Adding field 'Synthesis.product_smiles'
        db.add_column(u'app_synthesis', 'product_smiles',
                      self.gf('django.db.models.fields.CharField')(default='CCCCC', max_length=100),
                      keep_default=False)

        # Adding field 'Synthesis.correct_answer'
        db.add_column(u'app_synthesis', 'correct_answer',
                      self.gf('django.db.models.fields.CharField')(default='Magic', max_length=100),
                      keep_default=False)


    def backwards(self, orm):
        # Deleting model 'SingleStepHardProblem'
        db.delete_table(u'app_singlestephardproblem')

        # Adding field 'SingleStepProblem.productSmiles'
        db.add_column(u'app_singlestepproblem', 'productSmiles',
                      self.gf('django.db.models.fields.CharField')(default='CCCCC', max_length=100),
                      keep_default=False)

        # Adding field 'SingleStepProblem.reactantSmiles'
        db.add_column(u'app_singlestepproblem', 'reactantSmiles',
                      self.gf('django.db.models.fields.CharField')(default='CCC=CC', max_length=100),
                      keep_default=False)

        # Adding field 'SingleStepProblem.correctAnswer'
        db.add_column(u'app_singlestepproblem', 'correctAnswer',
                      self.gf('django.db.models.fields.CharField')(default='Magic', max_length=100),
                      keep_default=False)

        # Deleting field 'SingleStepProblem.reactant_smiles'
        db.delete_column(u'app_singlestepproblem', 'reactant_smiles')

        # Deleting field 'SingleStepProblem.product_smiles'
        db.delete_column(u'app_singlestepproblem', 'product_smiles')

        # Deleting field 'SingleStepProblem.correct_answer'
        db.delete_column(u'app_singlestepproblem', 'correct_answer')

        # Adding field 'PredictProductsProblem.reactantSmiles'
        db.add_column(u'app_predictproductsproblem', 'reactantSmiles',
                      self.gf('django.db.models.fields.CharField')(default='CCC=CC', max_length=100),
                      keep_default=False)

        # Adding field 'PredictProductsProblem.correctAnswer'
        db.add_column(u'app_predictproductsproblem', 'correctAnswer',
                      self.gf('django.db.models.fields.CharField')(default='Magic', max_length=100),
                      keep_default=False)

        # Deleting field 'PredictProductsProblem.reactant_smiles'
        db.delete_column(u'app_predictproductsproblem', 'reactant_smiles')

        # Deleting field 'PredictProductsProblem.correct_answer'
        db.delete_column(u'app_predictproductsproblem', 'correct_answer')

        # Deleting field 'PredictProductsProblem.incorrect_answer1'
        db.delete_column(u'app_predictproductsproblem', 'incorrect_answer1')

        # Deleting field 'PredictProductsProblem.incorrect_answer2'
        db.delete_column(u'app_predictproductsproblem', 'incorrect_answer2')

        # Deleting field 'PredictProductsProblem.incorrect_answer3'
        db.delete_column(u'app_predictproductsproblem', 'incorrect_answer3')

        # Adding field 'Synthesis.productSmiles'
        db.add_column(u'app_synthesis', 'productSmiles',
                      self.gf('django.db.models.fields.CharField')(default='BrCCl', max_length=100),
                      keep_default=False)

        # Adding field 'Synthesis.reactantSmiles'
        db.add_column(u'app_synthesis', 'reactantSmiles',
                      self.gf('django.db.models.fields.CharField')(default='CCCC=CC', max_length=100),
                      keep_default=False)

        # Adding field 'Synthesis.correctAnswer'
        db.add_column(u'app_synthesis', 'correctAnswer',
                      self.gf('django.db.models.fields.CharField')(default='H2 cat PdC', max_length=100),
                      keep_default=False)

        # Deleting field 'Synthesis.reactant_smiles'
        db.delete_column(u'app_synthesis', 'reactant_smiles')

        # Deleting field 'Synthesis.product_smiles'
        db.delete_column(u'app_synthesis', 'product_smiles')

        # Deleting field 'Synthesis.correct_answer'
        db.delete_column(u'app_synthesis', 'correct_answer')


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
            'correct_answer': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'incorrect_answer1': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'incorrect_answer2': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'incorrect_answer3': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'reactant_smiles': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'reagents': ('django.db.models.fields.CharField', [], {'max_length': '100'})
        },
        u'app.singlestephardproblem': {
            'Meta': {'object_name': 'SingleStepHardProblem'},
            'answer': ('django.db.models.fields.related.ForeignKey', [], {'to': u"orm['api.ReagentSet']"}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'product_smiles': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'reactant_smiles': ('django.db.models.fields.CharField', [], {'max_length': '100'})
        },
        u'app.singlestepproblem': {
            'Meta': {'object_name': 'SingleStepProblem'},
            'correct_answer': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'product_smiles': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'reactant_smiles': ('django.db.models.fields.CharField', [], {'max_length': '100'})
        },
        u'app.synthesis': {
            'Meta': {'object_name': 'Synthesis'},
            'correct_answer': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'id': ('django.db.models.fields.AutoField', [], {'primary_key': 'True'}),
            'product_smiles': ('django.db.models.fields.CharField', [], {'max_length': '100'}),
            'reactant_smiles': ('django.db.models.fields.CharField', [], {'max_length': '100'})
        }
    }

    complete_apps = ['app']