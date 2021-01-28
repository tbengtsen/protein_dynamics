# sys.path.insert(0, '/z/home/tebn/machine_learning_programs/gvp/src')
# from models import *
# # from 
# import util
# from datasets import parse_batch as parse_batch_gvp

# model = CPDModel(node_features=(8, 100), edge_features=(1,32), hidden_dim=(16,100))
# optimizer = tf.keras.optimizers.Adam()
# pretrained_gvp = '/z/home/tebn/machine_learning_programs/gvp/models/cath_pretrained.index'
# util.load_checkpoint(model, optimizer,pretrained_gvp ) 
# # predictor = model.sample


def _loss(S, log_probs, mask, num_letters=20):
    """ Negative log probabilities.
    from ingrahams test_rocklins mutations"""
    criterion = torch.nn.NLLLoss(reduction='none')
    loss = criterion(
        log_probs.contiguous().view(-1,num_letters),
        S.contiguous().view(-1)
    ).view(S.size())

    return loss

def tf_loss(S, logits, mask):
    loss_metric = tf.keras.metrics.SparseCategoricalCrossentropy(from_logits=True)
    loss_metric.update_state(S, logits, sample_weight=mask)
    loss  = loss_metric.result()
    
    return loss
    
def gvp_predictor(coords, seq, model):
    '''based on gvp predicts the probability of a given sequence 
    based on the log likelihood
    '''
    # convert to format that fits gvp
    batch = [{'seq': seq, 'coords': coords }]
    X, S, mask = parse_batch_gvp(batch) 
    
    # predict freq for each aa at each res position
    logits = model(X,S,mask,train=False)
    print(type(logits))
    
    loss =  tf_loss(S, logits, mask)
    print("loss",loss)
    # get approx to gibbs free energy -> -log(pi)=deltaGi => 
#     neglogP = torch.sum(loss * mask, dim=1) / torch.sum(mask, dim=1)
#     neglogP = neglogp.cpu().data.numpy().tolist()[0]
    
#     return neglogP

### LOOPS

# loop over all proteins
for protein in ssl.proteins[2:3]:
    print('seq len', len(ssl[protein]['seq'] ))
    print(f'\nCalculating predictions for {protein}:\n---------------------------------------')
    batch = [{'seq': ssl[protein]['seq'], 'coords':  ssl[protein]['coords'] }]
    
    X, S, mask = parse_batch_gvp(batch) 
    
    output = model(X,S,mask)
#     print(output.shape)
    seq =  ssl[protein]['seq']
    coords = ssl[protein]['coords'] 
    neglogP = gvp_predictor(coords, seq, model)
    print(neglogP)
#     design = util.sample(model, X, mask, 1)
#     print(design)
    
  
